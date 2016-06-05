import sys, os
import numpy as np

from ccdproc import CCDData
from pyhrs.hrsprocess import *

from astropy import units as u

import argparse


if __name__=='__main__':


   parser = argparse.ArgumentParser()
   parser.add_argument("-i","--infile",help="Name of the raw file to be processed",type=str,required=True)
   parser.add_argument("-b","--bias",help="Path to the bias frame (e.g. [mr]/HBIAS.fits)",type=str,required=False)
   parser.add_argument("-o","--orderframe",help="Path to the order frame (e.g. [mr]/HORDER.fits)",type=str,required=True)
   parser.add_argument("-f","--flat",help="Path to the flat field frame (e.g. [mr]/HFLAT.fits)",type=str,required=True)
   parser.add_argument("-n",'--noclean',help="The file will not be cleaned from cosmic rays", action='store_true')
   args=parser.parse_args()



   print args.bias
   if args.bias=='None' or args.bias==None: 
      bias=None
   else:
      bias = CCDData.read(args.bias, units=u.adu, ignore_missing_end=True)
   if os.path.basename(args.infile).startswith('H'):
        ccd = blue_process(args.infile, masterbias=bias)
   elif os.path.basename(args.infile).startswith('R'):
        ccd = red_process(args.infile, masterbias=bias)
   else:
        exit('Are you sure this is an HRS file?')


   if not args.noclean:
      from astroscrappy import detect_cosmics
      crmask, cleanarr = detect_cosmics(ccd.data, inmask=None, sigclip=4.5, sigfrac=0.3,
                      objlim=5.0, gain=1.0, readnoise=6.5,
                      satlevel=65536.0, pssl=0.0, niter=4,
                      sepmed=True, cleantype='meanmask', fsmode='median',
                      psfmodel='gauss', psffwhm=2.5, psfsize=7,
                      psfk=None, psfbeta=4.765, verbose=False) 
      ccd.data = cleanarr
   order_frame = CCDData.read(args.orderframe, unit=u.adu, ignore_missing_end=True)
   flat_frame = CCDData.read(args.flat, ignore_missing_end=True)
   ccd=flatfield_science(ccd, flat_frame, order_frame, median_filter_size=None, interp=True) 

   outfile = 'p'+args.infile
   ccd.write(outfile, clobber=True)

     
