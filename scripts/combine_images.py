import sys, os
import numpy as np
from astropy.io import fits



if __name__=='__main__':

   inlist = sys.argv[1]
   outfile = sys.argv[2]

   f=open(inlist,"r")
   
   #read fits files
   data_arr=[]
   for line in f.readlines():
      hdulist = fits.open(line.split("\n")[0])
      prihdr = hdulist[0].header
      ccd=hdulist[0].data
      data_arr.append(ccd)
      hdulist.close()
      
   #calculate mean
   out_arr=np.zeros_like(data_arr[0])
   i=0.
   for arr in data_arr:
      out_arr+=arr
      i+=1.
   out_arr=out_arr/i


#   import matplotlib.pyplot as plt
#   plt.imshow(out_arr, cmap='gray')
#   plt.colorbar()
#   plt.show()

   #write output
   hdu = fits.PrimaryHDU(out_arr)
   hdu.header=prihdr
   hdu.writeto(outfile, clobber=True)
