
import sys
from astropy.io import fits
import numpy as np
from pylab import *
from scipy import integrate

#figure()
#sys.argv[1] - sci fits; sys.argv[2] -cal fits ; sys.argv[3] - file cal std ; sys.argv[4] - optional filter ; sys.argv[5] - optional mag ; sys.argv[5] - optional zeropoint 

if __name__=='__main__':

   img_sci = sys.argv[1]
   hdu_sci = fits.open(img_sci)
   wave_sci = hdu_sci[1].data['Wavelength']
   flux_sci = hdu_sci[1].data['Flux']
   order_sci = hdu_sci[1].data['Order']

   i_sort=np.argsort(wave_sci)
   wave_sci=wave_sci[i_sort]
   flux_sci=flux_sci[i_sort]

   img_cal = sys.argv[2]
   hdu_cal = fits.open(img_cal)
   wave_cal = hdu_cal[1].data['Wavelength']
   flux_cal = hdu_cal[1].data['Flux']


   wave_min=np.min(wave_cal)
   wave_max=np.max(wave_cal)

   wave_std=[]
   flux_std=[]
   bin_std=[]

   #The 'f' files list wavelength ( A ), flux ( ergs/cm/cm/s/A * 10**16 )
   file_std=open(sys.argv[3])
   for line in file_std.readlines():
      wave_std.append(line.split()[0])
      flux_std.append(line.split()[1])
      bin_std.append(line.split()[3])
   wave_std=np.array(wave_std).astype(float)
   flux_std=np.array(flux_std).astype(float)
   bin_std=np.array(bin_std).astype(float)

   wave_ratio=[]
   ratio=[]

   for i in range(len(wave_std)):
      if (wave_std[i]-0.5*bin_std[i]>wave_min) and (wave_std[i]+0.5*bin_std[i]<wave_max):

         i_cal=np.where((wave_std[i]-0.5*bin_std[i]<wave_cal) & (wave_std[i]+0.5*bin_std[i]>wave_cal))
   #      int_cal=integrate.simps(flux_cal[i_cal],x=wave_cal[i_cal])
         int_cal=np.median(flux_cal[i_cal])
         wave_ratio.append(wave_std[i])
         ratio.append(int_cal/flux_std[i])

   ratio=np.array(ratio)
   wave_ratio=np.array(wave_ratio)
 

   check=True
   while check:
       check_input=True
       while check_input:
          check_input=False
          fit_order = raw_input("Order of the fit (default = 8): ")
          if fit_order=="":
             fit_order=6
          else:
             try:
                fit_order=int(fit_order)
             except ValueError:
                print "incorrect input"
                check_input=True

       check_input=True
       while check_input:
          check_input=False
          iterations = raw_input("Number of iterations (default = 5): ")
          if iterations=="":
             iterations=5
          else:
             try:
                iterations=int(iterations)
             except ValueError:
                print "incorrect input"
                check_input=True

       check_input=True
       while check_input:
          check_input=False
          lowrej = raw_input("Lower rejection limit (default = 3): ")
          if lowrej=="":
             lowrej=3
          else:
             try:
                lowrej=float(lowrej)
             except ValueError:
                print "incorrect input"
                check_input=True

       check_input=True
       while check_input:
          check_input=False
          uprej = raw_input("Upper rejection limit (default = 3): ")
          if uprej=="":
             uprej=3
          else:
             try:
                uprej=float(uprej)
             except ValueError:
                print "incorrect input"
                check_input=True
      
       for i in range(iterations):
          if i==0:
             coefficients=np.polyfit(wave_ratio,ratio, fit_order)
          else:

             if i_fit[0].size==0:
                print "error while rejecting point - try higher values of uprej and lowrej"
                i_fit=np.arange(ratio)

             coefficients=np.polyfit(wave_ratio[i_fit],ratio[i_fit], fit_order)

          func=np.poly1d(coefficients)
          fitted=np.polyval(func,wave_ratio)
          residuals=(ratio-fitted)

          mean_res=np.mean(residuals)
          std_res=np.std(residuals)
          i_fit=np.where(  (residuals<uprej*std_res)  & (residuals>-lowrej*std_res))
 

       fitted_sci=np.polyval(func,wave_sci)
       flux_s=flux_sci/fitted_sci
      
       figure()
       subplot(2,1, 1)
       scatter(wave_ratio,ratio)
       plot(wave_sci,fitted_sci)
       ylabel("Ratios")
       subplot(2, 1, 2)
       scale_sci=np.median(flux_sci)
       scale_s=np.median(flux_s)
       plot(wave_sci,flux_sci/scale_sci)
       plot(wave_sci,flux_s/scale_s)
       xlabel(r"Wavelength [$\AA$]")
       ylabel("Normalized spectrum ")
       show()


       check_input=True
       while check_input==True:
          repeat=raw_input("Repeat fitting? (default = no): ")
          if repeat=="y" or repeat=="yes" or repeat=="Y" or repeat=="YES":
             check_input=False
          elif repeat=="n" or repeat=="no" or repeat=="N" or repeat=="NO" or repeat=="":
             check_input=False
             check=False
       print
   


   

   c1 = fits.Column(name='Wavelength', format='D', array=wave_sci, unit='Angstroms')
   c2 = fits.Column(name='Flux', format='D', array=flux_s, unit='Counts')
   c3 = fits.Column(name='Order', format='I', array=order_sci)


   outfile="f"+sys.argv[1]
   tbhdu = fits.BinTableHDU.from_columns([c1,c2,c3])
   tbhdu.writeto(outfile, clobber=True)

