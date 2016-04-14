import sys, os
import numpy as np
from astropy.io import fits

from pyhrs import create_masterbias, create_masterflat, create_orderframe, wavelength_calibrate_arc, wavelength_calibrate_order
      

import argparse





if __name__=='__main__':

   parser = argparse.ArgumentParser()
   parser.add_argument("fits_list",help="Text file with list of HRS observations ",type=str)
 #  parser.add_argument("-a","--all",help="Fits all orders simultaneously",action="store_true")
   args=parser.parse_args()



   lr_Hflats=[]
   lr_Rflats=[]
   lr_Harcs=[]
   lr_Rarcs=[]
   mr_Hflats=[]
   mr_Rflats=[]
   mr_Harcs=[]
   mr_Rarcs=[]
   hr_Hflats=[]
   hr_Rflats=[]
   hr_Harcs=[]
   hr_Rarcs=[]
   Hbiases=[]
   Rbiases=[]

   f=open(args.fits_list,"r")

   for line in f.readlines():
      if not (line[0]=="H" or line[0]=="R"):
         print "Are you sure this is HRS observation? ",line
      else:
         fits_name=line.split("\n")[0]
         hdulist = fits.open(fits_name)
         mode=hdulist[0].header['OBSMODE']
         obj=hdulist[0].header['OBJECT']
         hdulist.close()
         

         if line[0]=="H":
            if mode=="LOW RESOLUTION":
               if obj=="Arc":
                  lr_Harcs.append(fits_name)
               elif obj=="Flat field":
                  lr_Hflats.append(fits_name)
            elif mode=="MEDIUM RESOLUTION":
               if obj=="Arc":
                  mr_Harcs.append(fits_name)
               elif obj=="Flat field":
                  mr_Hflats.append(fits_name)
            elif mode=="HIGH RESOLUTION":
               if obj=="Arc":
                  hr_Harcs.append(fits_name)
               elif obj=="Flat field":
                  hr_Hflats.append(fits_name)
            elif obj=="Bias":
               Hbiases.append(fits_name)
         elif line[0]=="R":
            if mode=="LOW RESOLUTION":
               if obj=="Arc":
                  lr_Rarcs.append(fits_name)
               elif obj=="Flat field":
                  lr_Rflats.append(fits_name)
            elif mode=="MEDIUM RESOLUTION":
               if obj=="Arc":
                  mr_Rarcs.append(fits_name)
               elif obj=="Flat field":
                  mr_Rflats.append(fits_name)
            elif mode=="HIGH RESOLUTION":
               if obj=="Arc":
                  hr_Rarcs.append(fits_name)
               elif obj=="Flat field":
                  hr_Rflats.append(fits_name)
            elif obj=="Bias":
               Rbiases.append(fits_name)
   f.close()

   if not len(lr_Harcs)==0:
      print "ard"
