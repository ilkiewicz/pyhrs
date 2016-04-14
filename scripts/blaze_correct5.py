
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


   of=np.ones_like(fit_wave[i_fit])*float(o)
   mw=np.ones_like(fit_wave[i_fit])*np.mean(fit_wave[i_fit])
   wf_arr=np.append(wf_arr,fit_wave[i_fit])
   mf_arr=np.append(mf_arr,mw)
   rf_arr=np.append(rf_arr,(fit_ratio[i_fit]/np.mean(fit_ratio[i_fit])))
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

