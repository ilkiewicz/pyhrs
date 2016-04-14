
import sys
from astropy.io import fits
import numpy as np
import glob
from pylab import *






figure()

img = sys.argv[1]
hdu = fits.open(img)
wave = hdu[1].data['Wavelength']
flux = hdu[1].data['Flux']
order = hdu[1].data['Order']

if (sys.argv[2]=='all') or (sys.argv[2]=='ALL') or (sys.argv[2]=='All'):
   for o in np.unique(order):
      o = int(o)
      mask = (order==o)
      flux[mask] = np.convolve(flux[mask], np.ones(10), 'same')
      plot(wave[mask], flux[mask])
else:
   for o in sys.argv[2:]:
      o = int(o)
      mask = (order==o)
      flux[mask] = np.convolve(flux[mask], np.ones(10), 'same')
      plot(wave[mask], flux[mask])

show()
