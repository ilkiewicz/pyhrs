
import sys
from astropy.io import fits
import numpy as np
from specutils.io import write_fits

#figure()
#sys.argv[1] - sci fits


img_sci = sys.argv[1]
hdu_sci = fits.open(img_sci)
wave_sci = hdu_sci[1].data['Wavelength']
flux_sci = hdu_sci[1].data['Flux']
order_sci = hdu_sci[1].data['Order']


from specutils.io import write_fits


if sys.argv[2]=="iraf":

   outfile=sys.argv[1].replace("spec","iraf")
   spec1d = Spectrum1D.from_array(wave_sci, flux_sci)
   write_fits.write(spec1d, outfile)
elif sys.argv[2]=="ascii" or sys.argv[2]=="ASCII" or sys.argv[2]=="dat":
   outfile=sys.argv[1].replace("fits","dat")
   f=open(outfile,"w")
   for i in range(len(wave_sci)):
      f.write(str(wave_sci[i])+" "+str(flux_sci[i])+" "+str(order_sci[i])+"\n")
   f.close()
