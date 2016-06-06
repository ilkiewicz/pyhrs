
import sys
from astropy.io import fits
import numpy as np
from pylab import *
from scipy import integrate
import math
import argparse





def merge_orders(wave,flux,order,keep=False):
   """Scales echelle orders using overlapping parts and merges the result. By default the result wavelength table will contain only one measurement per wavelength (i.e. in the overlapping part only one measurement will be kept). 

   Parameters
   ----------
   wave: array of floats
       Array containing wavelength
   flux: array of floats
       Array containing measured counts/flux
   order: array of floats
       Array containing numbers of orders of corresponding wave/flux
   keep: bool
       If true the overlapping part will not be deleted

   Returns
   -------
   w_arr: array of floats
       Array containing wavelength
   f_arr: array of floats
       Array containing scaled counts/flux
   o_arr: array of floats
       Array containing order of corresponding wave/flux

   """
   w_arr=np.array([])
   f_arr=np.array([])
   o_arr=np.array([])


   list_orders=np.unique(order)

   for o in list_orders:
      o=int(o)
      if o==int(list_orders[0]):    #the first order in fits
         o_l1=np.where(order==o)  
         o_l2=np.where(order==(o+1))
         w_tmp1=wave[o_l1]
         f_tmp1=flux[o_l1] 
         o_tmp1=order[o_l1] 
         w_tmp2=wave[o_l2]
         w1=np.mean(w_tmp1)  
         w2=np.mean(w_tmp2)  
         if keep:
            w_arr=np.append(w_arr,w_tmp1) 
            f_arr=np.append(f_arr,f_tmp1) 
            o_arr=np.append(o_arr,o_tmp1)
         else:
            for i in range(len(w_tmp1)):
               if math.fabs(w_tmp1[i]-w1)< math.fabs(w_tmp1[i]-w2):
                  w_arr=np.append(w_arr,[w_tmp1[i]]) 
                  f_arr=np.append(f_arr,[f_tmp1[i]]) 
                  o_arr=np.append(o_arr,[o_tmp1[i]])

      elif o==int(list_orders[len(list_orders)-1]):    #the last order in fits
         o_l0=np.where(order==o-1)  
         o_l1=np.where(order==o)  
         w_tmp0=wave[o_l0]
         w_tmp1=wave[o_l1]
         f_tmp0=flux[o_l0]
         f_tmp1=flux[o_l1] 
         o_tmp1=order[o_l1] 
         w0=np.mean(w_tmp0) 
         w1=np.mean(w_tmp1)  
         
         
         min_overlap=np.min(w_arr)+0.15*(np.max(w_tmp1)-np.min(w_arr))
         max_overlap=np.max(w_tmp1)-0.15*(np.max(w_tmp1)-np.min(w_arr))


         i_1=np.where((w_tmp1>min_overlap) & (w_tmp1<max_overlap))
         i_a=np.where((w_arr>min_overlap) & (w_arr<max_overlap))
         f_tmp1=f_tmp1*np.median(f_arr[i_a])/np.median(f_tmp1[i_1])
         if keep:
            w_arr=np.append(w_arr,w_tmp1) 
            f_arr=np.append(f_arr,f_tmp1) 
            o_arr=np.append(o_arr,o_tmp1)
         else:
            for i in range(len(w_tmp1)):
               if (math.fabs(w_tmp1[i]-w1)< math.fabs(w_tmp1[i]-w0)):
                  w_arr=np.append(w_arr,[w_tmp1[i]]) 
                  f_arr=np.append(f_arr,[f_tmp1[i]])  
                  o_arr=np.append(o_arr,[o_tmp1[i]])

      else:    #every other order in fits
         o_l0=np.where(order==o-1)  
         o_l1=np.where(order==o)  
         o_l2=np.where(order==(o+1))
         w_tmp0=wave[o_l0]
         w_tmp1=wave[o_l1]
         f_tmp0=flux[o_l0]
         f_tmp1=flux[o_l1] 
         o_tmp1=order[o_l1] 
         w_tmp2=wave[o_l2]
         w0=np.mean(w_tmp0) 
         w1=np.mean(w_tmp1)  
         w2=np.mean(w_tmp2) 
         min_overlap=np.min(w_arr)+0.2*(np.max(w_tmp1)-np.min(w_arr))
         max_overlap=np.max(w_tmp1)-0.05*(np.max(w_tmp1)-np.min(w_arr))

         i_1=np.where((w_tmp1>min_overlap) & (w_tmp1<max_overlap))
         i_a=np.where((w_arr>min_overlap) & (w_arr<max_overlap))
         f_tmp1=f_tmp1*np.median(f_arr[i_a])/np.median(f_tmp1[i_1])


         if keep:
            w_arr=np.append(w_arr,w_tmp1) 
            f_arr=np.append(f_arr,f_tmp1) 
            o_arr=np.append(o_arr,o_tmp1)
         else:
            for i in range(len(w_tmp1)):
               if (math.fabs(w_tmp1[i]-w1)< math.fabs(w_tmp1[i]-w2)) and (math.fabs(w_tmp1[i]-w1)< math.fabs(w_tmp1[i]-w0)):
                  w_arr=np.append(w_arr,[w_tmp1[i]]) 
                  f_arr=np.append(f_arr,[f_tmp1[i]]) 
                  o_arr=np.append(o_arr,[o_tmp1[i]]) 

   return w_arr,f_arr,o_arr




if __name__=='__main__':

   parser = argparse.ArgumentParser()
   parser.add_argument("-i","--infile",help="Name of the file with spectrum to be merged",type=str,required=True)
   parser.add_argument("-k","--keep",help="The overlapping regions will not be deleted", action='store_true')
   args=parser.parse_args()

   img = args.infile
   hdu = fits.open(img)
   wave = hdu[1].data['Wavelength']
   flux = hdu[1].data['Flux']
   order = hdu[1].data['Order']



   w_arr,f_arr,o_arr=merge_orders(wave,flux,order,keep=args.keep)



   c1 = fits.Column(name='Wavelength', format='D', array=w_arr, unit='Angstroms')
   c2 = fits.Column(name='Flux', format='D', array=f_arr, unit='Counts')
   c3 = fits.Column(name='Order', format='I', array=o_arr)

   outfile="m"+args.infile
   tbhdu = fits.BinTableHDU.from_columns([c1,c2,c3])
   tbhdu.writeto(outfile, clobber=True)

