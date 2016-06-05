
import sys
from astropy.io import fits
import numpy as np
from pylab import *

from astropy import modeling as amod
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
m_arr=np.array([])
f_arr=np.array([])
o_arr=np.array([])

wf_arr=np.array([])
rf_arr=np.array([])
of_arr=np.array([])
mf_arr=np.array([])

for o in np.unique(order_sci):
   o = int(o)

   order_i=np.where(order_sci==o )
   order_j=np.where(order_flat==o )




   #selecting points to fit - deleting outliers

   #fitting low order function to flat spectrum
   fit_order=1
   coefficients=np.polyfit(wave_flat[order_j], flux_flat[order_j], fit_order)
   func=np.poly1d(coefficients)

   fitted=np.polyval(func,wave_flat)
#   plot(wave_flat[order_j],fitted)
#   scatter(wave_flat[order_j], flux_flat[order_j])
#   show()

   #selecting outliers

   residuals=np.abs(flux_flat-fitted)
   mean_res=np.mean(residuals[order_j])
   std_res=np.std(residuals[order_j])
 

#   tmp=np.where((residuals < mean_res+4.*std_res) & (wave_flat<np.max(wave_sci[order_i])) & (wave_flat>np.min(wave_sci[order_i]))& (order_flat==o))
#   scatter(wave_flat[order_j],residuals[order_j])
#   scatter(wave_flat[tmp],residuals[tmp],c='r')
#   show()

 
   bins=50.

   wave_min=np.min(wave_sci[order_i])
   wave_max=np.max(wave_sci[order_i])
   step=(wave_max-wave_min)/bins
   fit_wave=[]
   fit_ratio=[]
   for wv in np.arange(wave_min,wave_max,step):
      i_flat=np.where((wave_flat >= wv) & (wave_flat <= (wv+step)) & (order_flat==o) & (residuals < mean_res+3.*std_res) ) #check which elements of array wave_flat are within spectral range "wv" and "wv+step"
      if len(i_flat[0])==0:
         i_flat=np.where((wave_flat >= wv) & (wave_flat <= (wv+step)) & (order_flat==o) & (residuals < 1e9) )
         
#      int_flux=integrate.simps(flux_flat[i_flat],x=wave_flat[i_flat]) 
      int_flux=np.median(flux_flat[i_flat])
      i_cal=np.where((wave_cal >= wv) & (wave_cal <= (wv+step)) )  #check which elements of array wave_flat are within spectral range "wv" and "wv+step"
#      int_cal=integrate.simps(flux_cal[i_cal],x=wave_cal[i_cal])
      int_cal=np.median(flux_cal[i_cal])
   
      fit_wave.append(wv+0.5*step)
      fit_ratio.append(int_flux/int_cal)

   fit_wave=np.array(fit_wave)
   fit_ratio=np.array(fit_ratio)  


   #first estimate fit to ratios
   fit_order=2
   coefficients=np.polyfit(fit_wave, fit_ratio, fit_order)
   func=np.poly1d(coefficients)
   fitted=np.polyval(func,fit_wave) 


   residuals=np.abs(fitted-fit_ratio)
   mean_fit=np.mean(residuals)
   std_fit=np.std(residuals)


   mean=np.mean(fit_ratio)
   std=np.std(fit_ratio)
   residuals_m=np.abs(fit_ratio-mean)#residuals from a mean
 
   i_fit=np.where((residuals < mean_fit+3.*std_fit) & (residuals_m < mean+3.*std)) 

   if len(i_fit[0])==0:
      i_fit=np.where(residuals < 1.e9) 

   #finatl fit to ratios
#   fit_order=6
#   coefficients=np.polyfit(fit_wave[i_fit], fit_ratio[i_fit], fit_order)
#   func=np.poly1d(coefficients)
#   fitted=np.polyval(func,wave_sci)

   of=np.ones_like(fit_wave[i_fit])*float(o)
   mw=np.ones_like(fit_wave[i_fit])*np.mean(fit_wave[i_fit])
   wf_arr=np.append(wf_arr,fit_wave[i_fit])
   mf_arr=np.append(mf_arr,mw)
   rf_arr=np.append(rf_arr,(fit_ratio[i_fit]/np.mean(fit_ratio[i_fit])))
#   rf_arr=np.append(rf_arr,(fit_ratio[i_fit]))
   of_arr=np.append(of_arr,of)

   w_arr=np.append(w_arr,wave_sci[order_i])
   m=np.ones_like(wave_sci[order_i])*np.mean(wave_sci[order_i])
   m_arr=np.append(m_arr,m)
   f_arr=np.append(f_arr,flux_sci[order_i])
   o_arr=np.append(o_arr,order_sci[order_i])

p_init = amod.models.Polynomial2D(degree=20)
fit_p = amod.fitting.LevMarLSQFitter()
p = fit_p(p_init,wf_arr-mf_arr, of_arr, rf_arr)

figure()
subplot(1, 3, 1)
scatter(wf_arr-mf_arr, of_arr, c=rf_arr,edgecolor="None")
xlabel("X")
ylabel("ORDER")
colorbar(orientation='horizontal')
title("Ratios")
subplot(1, 3, 2)
scatter(wf_arr-mf_arr, of_arr, c=p(wf_arr-mf_arr, of_arr),edgecolor="None")
xlabel("X")
colorbar(orientation='horizontal')
title("2D Fit")
subplot(1, 3, 3)
scatter(wf_arr-mf_arr, of_arr, c=(rf_arr-p(wf_arr-mf_arr, of_arr)),edgecolor="None")
xlabel("X")
colorbar(orientation='horizontal')
title("Residual")
show()


c1 = fits.Column(name='Wavelength', format='D', array=w_arr, unit='Angstroms')
c2 = fits.Column(name='Flux', format='D', array=(f_arr/p(w_arr-m_arr, o_arr)), unit='Counts')
c3 = fits.Column(name='Order', format='I', array=o_arr)

outfile="b"+sys.argv[1]
tbhdu = fits.BinTableHDU.from_columns([c1,c2,c3])
tbhdu.writeto(outfile, clobber=True)
