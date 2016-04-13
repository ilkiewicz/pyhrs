
import sys
from astropy.io import fits
import numpy as np
import pylab as pl
import argparse

def plot_fit(warr,farr,fitted,i_fit):

   #making list of outliers
   i_outliers=range(len(farr))
   i_all=np.where(farr>0.)
   mask1=np.in1d(i_outliers, i_fit)
   mask2=np.in1d(i_outliers, i_all)
   i_outliers=np.where(~mask1 & mask2)

   figure = pl.figure()
   ax1 = figure.add_subplot(211)
   ax1.plot(warr[i_all],farr[i_all],c="green") 
   ax1.scatter(warr[i_outliers],farr[i_outliers],c="red",edgecolor="None") 
      
   ax1.axes.get_xaxis().set_visible(False)
   ax1.plot(warr[i_all],fitted[i_all],c="blue")

   ax2 = figure.add_subplot(212, sharex=ax1)
   ax2.hold(False)
   ax2.plot(warr[i_all],farr[i_all]/fitted[i_all],c="blue")
   ax2.set_xlabel(r"Wavelength [$\AA$]")

   pl.tight_layout()
   pl.show()


def fit_spec(warr,farr,iterations,fit_order,lowrej,uprej,interactive=False):
   for i in range(iterations):          #iterative fitting
      if i==0:                          #first iteration - there are no residuals available
         coefficients=np.polyfit(warr, farr, fit_order)
         i_fit=i_fit=np.where(farr > 0.)
        
      else:                             #fitting only to good points (not outliers)
         if i_fit[0].size==0:           #if all point are outliers print error
            print "Error while rejecting points - try higher values of upper rejection and lower rejection"
            i_fit=np.where(farr > 0.)
         coefficients=np.polyfit(warr[i_fit], farr[i_fit], fit_order)
               

      func=np.poly1d(coefficients)
      fitted=np.polyval(func,warr)
      residuals=(farr-fitted)


      #rejecting outliers and bad regions
      mean_res=np.mean(residuals[i_fit])
      std_res=np.std(residuals[i_fit])

      i_fit=np.where( (farr > 0.) & (residuals<uprej*std_res)  & (residuals>-lowrej*std_res))
      
   if interactive:
      plot_fit(warr,farr,fitted,i_fit)

   print "Mean residual of fit: ",  mean_res     
   print "Standard deviation of fit: ",  std_res
   return coefficients,i_fit,fitted



def get_fit_params(fit_order_d=3,iterations_d=15,lowrej_d=1.5,uprej_d=3):
  check_input=True
  while check_input:
     check_input=False
     fit_order = raw_input("Order of the fit (default = %i): "%int(fit_order_d))
     if fit_order=="":
        fit_order=int(fit_order_d)
     else:
        try:
           fit_order=int(fit_order)
        except ValueError:
           print "incorrect input"
           check_input=True

  check_input=True
  while check_input:
     check_input=False
     iterations = raw_input("Number of iterations (default = %i): "%int(iterations_d))
     if iterations=="":
        iterations=int(iterations_d)
     else:
        try:
           iterations=int(iterations)
        except ValueError:
           print "incorrect input"
           check_input=True

  check_input=True
  while check_input:
     check_input=False
     lowrej = raw_input("Lower rejection limit (default = %.1f): "%float(lowrej_d))
     if lowrej=="":
        lowrej=float(lowrej_d)
     else:
        try:
           lowrej=float(lowrej)
        except ValueError:
           print "incorrect input"
           check_input=True

  check_input=True
  while check_input:
     check_input=False
     uprej = raw_input("Upper rejection limit (default = %.1f): "%float(uprej_d))
     if uprej=="":
        uprej=float(uprej_d)
     else:
        try:
           uprej=float(uprej)
        except ValueError:
           print "incorrect input"
           check_input=True

  return fit_order,iterations,lowrej,uprej


if __name__=='__main__':

   parser = argparse.ArgumentParser()
   parser.add_argument("spectrum_fits",help="Fits file with an extracted HRS spectrum",type=str)
   parser.add_argument('-f','--fit_pars', nargs='+', help='List of fitting parameters: fit_order, iterations, lowrej, uprej', required=False)
   parser.add_argument("-a","--all",help="Fits all orders simultaneously",action="store_true")
   parser.add_argument("-i","--interactive",help="Interactive fitting",action="store_true")
   args=parser.parse_args()

   if (args.interactive) or (not args.fit_pars==None):
      check_input=True
   else:
      print "Fitting parameters should be provided ('-f') or script should be run in interactive mode ('-i')"
      check_input=False
      
   if (not args.fit_pars==None):
      if len(args.fit_pars)==4:
         fit_order=int(args.fit_pars[0])
         iterations=int(args.fit_pars[1])
         lowrej=float(args.fit_pars[2])
         uprej=float(args.fit_pars[3])
      else:
         print "Error: list of fitting parameters should contain [fitting_order, interation, lower_rejection_lim, upper_rej_lim]"
         check_input=False

   if check_input:
      img_sci = args.spectrum_fits
      hdu_sci = fits.open(img_sci)
      wave_sci = hdu_sci[1].data['Wavelength']
      flux_sci = hdu_sci[1].data['Flux']
      order_sci = hdu_sci[1].data['Order']

      #sorting arrays to have 
      i_sort=np.argsort(wave_sci)
      wave_sci=wave_sci[i_sort]
      flux_sci=flux_sci[i_sort]
      order_sci=order_sci[i_sort]


      w_arr=np.array([])
      f_arr=np.array([])
      o_arr=np.array([])

      if not args.all:
         for o in np.unique(order_sci):
            o = int(o)
            order_i=np.where(order_sci==o )

            print 
            print "Normalizing order ",o

            check=True
            while check:
               check=False

               if args.interactive:
                  if not args.fit_pars==None:
                     fit_order,iterations,lowrej,uprej=get_fit_params(fit_order,iterations,lowrej,uprej)
                  else:
                     fit_order,iterations,lowrej,uprej=get_fit_params()          

               warr_o=wave_sci[order_i]
               farr_o=flux_sci[order_i]

               if args.interactive:
                  coefficients,i_fit,fitted=fit_spec(warr_o,farr_o,iterations,fit_order,lowrej,uprej,interactive=True)
               else:
                  coefficients,i_fit,fitted=fit_spec(warr_o,farr_o,iterations,fit_order,lowrej,uprej,interactive=False)
              
               if args.interactive:
                  check_input=True
                  while check_input==True:
                     repeat=raw_input("Repeat fitting? (default = no): ")
                     if repeat=="y" or repeat=="yes" or repeat=="Y" or repeat=="YES":
                        check_input=False
                     elif repeat=="n" or repeat=="no" or repeat=="N" or repeat=="NO" or repeat=="":
                        check_input=False
                        check=False
                  print
            w_arr=np.append(w_arr,warr_o)
            f_arr=np.append(f_arr,farr_o/fitted)
            o_arr=np.append(o_arr,order_sci[order_i])

      else:

         print 
         print "Normalizing whole spectrum "

         check=True
         while check:
            check=False

            if args.interactive:
               if not args.fit_pars==None:
                  fit_order,iterations,lowrej,uprej=get_fit_params(fit_order,iterations,lowrej,uprej)
               else:
                  fit_order,iterations,lowrej,uprej=get_fit_params()          

            if args.interactive:
               coefficients,i_fit,fitted=fit_spec(wave_sci,flux_sci,iterations,fit_order,lowrej,uprej,interactive=True)
            else:
               coefficients,i_fit,fitted=fit_spec(wave_sci,flux_sci,iterations,fit_order,lowrej,uprej,interactive=False)
            
            if args.interactive:
               check_input=True
               while check_input==True:
                  repeat=raw_input("Repeat fitting? (default = no): ")
                  if repeat=="y" or repeat=="yes" or repeat=="Y" or repeat=="YES":
                     check_input=False
                  elif repeat=="n" or repeat=="no" or repeat=="N" or repeat=="NO" or repeat=="":
                     check_input=False
                     check=False
               print
         w_arr=wave_sci
         f_arr=flux_sci/fitted
         o_arr=order_sci

      outfile="n"+args.spectrum_fits

      #Wavelength and Order columns need to be overwritten since the data was sorted
      hdu_sci[1].data['Wavelength']=w_arr
      hdu_sci[1].data['Flux']=f_arr
      hdu_sci[1].data['Order']=o_arr

      hdu_sci.writeto(outfile, clobber=True)

