
import sys
from astropy.io import fits
import numpy as np
from pylab import *
#from scipy import integrate

#figure()
#sys.argv[1] - sci fits; sys.argv[2] - flat file; sys.argv[3] - flat calib spec


img_sci = sys.argv[1]
hdu_sci = fits.open(img_sci)
wave_sci = hdu_sci[1].data['Wavelength']
flux_sci = hdu_sci[1].data['Flux']
order_sci = hdu_sci[1].data['Order']



img_flat = sys.argv[2]
hdu_flat = fits.open(img_flat)
wave_flat = hdu_flat[1].data['Wavelength']
flux_flat = hdu_flat[1].data['Flux']
order_flat = hdu_flat[1].data['Order']



spec_cal=open(sys.argv[3])

wave_cal=[]
flux_cal=[]
for line in spec_cal.readlines():
   wave_cal.append(line.split()[0])
   flux_cal.append(line.split()[1])
spec_cal.close()

wave_cal=np.array(wave_cal).astype(float)
flux_cal=np.array(flux_cal).astype(float)

w_arr=np.array([])
f_arr=np.array([])
o_arr=np.array([])

for o in np.unique(order_sci):
   o = int(o)

   order_i=np.where(order_sci==o )
   order_j=np.where(order_flat==o )



   #fitting low order function to a Flat spectrum for rejecting outliers
   fit_order=1
   coefficients=np.polyfit(wave_flat[order_j], flux_flat[order_j], fit_order)
   func=np.poly1d(coefficients)

   fitted_flat=np.polyval(func,wave_flat)
#   plot(wave_flat[order_j],fitted)
#   scatter(wave_flat[order_j], flux_flat[order_j])
#   show()

   fit_order=20
   iterations=2
   for i in range(iterations):
      #selecting outliers
      residuals=np.abs(flux_flat-fitted_flat)
      mean_res=np.mean(residuals[order_j])
      std_res=np.std(residuals[order_j])
      i_flat=np.where( (order_flat==o) & (residuals < mean_res+3.*std_res) ) #check which elements of array wave_flat are within spectral range "wv" and "wv+step"
      
      #fitting
      coefficients=np.polyfit(wave_flat[i_flat], flux_flat[i_flat], fit_order)
      func_flat=np.poly1d(coefficients)
      fitted_flat=np.polyval(func_flat,wave_flat)
 #  plt.plot(wave_flat[i_flat],fitted_flat[i_flat])
 #  plt.plot(wave_flat[i_flat],flux_flat[i_flat])
  # plt.show()


   
   #fitting function to a calibration spectrum

   #selecting part of the calibration spectrum
   wave_min=np.min(wave_sci[order_i])
   wave_max=np.max(wave_sci[order_i])
   i_cal=np.where( (wave_cal>wave_min) & (wave_cal<wave_max) )

   #fitting low order function to a calibration spectrum for rejecting outliers
   fit_order=8
   coefficients=np.polyfit(wave_cal[i_cal], flux_cal[i_cal], fit_order)
   func_cal=np.poly1d(coefficients)
 #  fitted_cal=np.polyval(func_cal,wave_cal)





   f_cal=np.polyval(func_cal,wave_sci[order_i])
   f_flat=np.polyval(func_flat,wave_sci[order_i])
   ratios=f_flat/f_cal

#   plot(wave_sci[order_i],flux_sci[order_i]/ratios)
#   show()



   plt.figure()
   plt.subplot(2, 1, 1)
   plt.plot(wave_flat[i_flat],fitted_flat[i_flat])
   plt.plot(wave_flat[i_flat],flux_flat[i_flat])
   plt.subplot(2, 1, 2)
   plot(wave_sci[order_i],flux_sci[order_i]/ratios)
   plt.show()

   w_arr=np.append(w_arr,wave_sci[order_i])
   f_arr=np.append(f_arr,flux_sci[order_i]/ratios)
   o_arr=np.append(o_arr,order_sci[order_i])

 #  plot(wave_sci[order_i],flux_sci[order_i]/fitted[order_i])
#show()

c1 = fits.Column(name='Wavelength', format='D', array=w_arr, unit='Angstroms')
c2 = fits.Column(name='Flux', format='D', array=f_arr, unit='Counts')
c3 = fits.Column(name='Order', format='I', array=o_arr)

outfile="b"+sys.argv[1]
tbhdu = fits.BinTableHDU.from_columns([c1,c2,c3])
tbhdu.writeto(outfile, clobber=True)

