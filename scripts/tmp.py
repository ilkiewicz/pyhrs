
import sys
from astropy.io import fits
import numpy as np
from pylab import *
from scipy import integrate
import math

#figure()
#sys.argv[1] - sci fits

img = sys.argv[1]
hdu = fits.open(img)
wave = hdu[1].data['Wavelength']
flux = hdu[1].data['Flux']
order = hdu[1].data['Order']



w_arr=np.array([])
f_arr=np.array([])
o_arr=np.array([])

#print np.unique(order)

"""
#making list of all orders except the last
list_orders=np.unique(order)
mask = np.ones(list_orders.shape,dtype=bool)
mask[len(mask)-1]=0
list_orders=list_orders[mask]
"""
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
      if "keep" in sys.argv:
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
      f_tmp1=f_tmp1
      if "keep" in sys.argv:
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
      f_tmp1=f_tmp1
     # plot(w_tmp0,f_tmp0)
     # plot(w_tmp1,f_tmp1)
    #  plot(w_tmp1,f_tmp1*np.median(f_tmp0[i_0])/np.median(f_tmp1[i_1]))
   #   plot(w_tmp0,w_tmp0*0.+np.median(f_tmp0[i_0]))
  #    plot(w_tmp1,w_tmp1*0.+np.median(f_tmp1[i_1]))
      #plot(w_tmp0[i_0],f_tmp0[i_0])
      #plot(w_tmp1[i_1],f_tmp1[i_1])
 #     show()
#      clf()
#      print np.median(f_tmp0[i_0]), np.median(f_tmp1[i_1])
      if "keep" in sys.argv:
         w_arr=np.append(w_arr,w_tmp1) 
         f_arr=np.append(f_arr,f_tmp1) 
         o_arr=np.append(o_arr,o_tmp1)
      else:
         for i in range(len(w_tmp1)):
            if (math.fabs(w_tmp1[i]-w1)< math.fabs(w_tmp1[i]-w2)) and (math.fabs(w_tmp1[i]-w1)< math.fabs(w_tmp1[i]-w0)):
               w_arr=np.append(w_arr,[w_tmp1[i]]) 
               f_arr=np.append(f_arr,[f_tmp1[i]]) 
               o_arr=np.append(o_arr,[o_tmp1[i]]) 


#   w_arr=np.append(w_arr,wave_sci[order_i])
#   f_arr=np.append(f_arr,flux_sci[order_i]/fitted[order_i])
#   o_arr=np.append(o_arr,order_sci[order_i])

#print w_arr.argsort()

c1 = fits.Column(name='Wavelength', format='D', array=w_arr, unit='Angstroms')
c2 = fits.Column(name='Flux', format='D', array=f_arr, unit='Counts')
c3 = fits.Column(name='Order', format='I', array=o_arr)

outfile="m"+sys.argv[1]
tbhdu = fits.BinTableHDU.from_columns([c1,c2,c3])
tbhdu.writeto(outfile, clobber=True)

