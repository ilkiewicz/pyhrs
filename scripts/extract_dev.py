import os
import sys
import numpy as np
import pickle

from ccdproc import CCDData

import specutils
from astropy import units as u

from astropy import modeling as mod
from astropy.io import fits

import pylab as pl

import specreduce

from specreduce.interidentify import InterIdentify
from specreduce import spectools as st
from specreduce import WavelengthSolution


from pyhrs import mode_setup_information
from pyhrs import zeropoint_shift
from pyhrs import HRSOrder, HRSModel

from numpy.polynomial import polynomial

def write_spdict(outfile, sp_dict):
    
    o_arr = None
    w_arr = None
    f_arr = None

    for k in sp_dict.keys():
        w,f = sp_dict[k]
        if w_arr is None:
            w_arr = 1.0*w
            f_arr = 1.0*f
            o_arr = k*np.ones_like(w, dtype=int)
        else:
            w_arr = np.concatenate((w_arr, w))
            f_arr = np.concatenate((f_arr, f))
            o_arr = np.concatenate((o_arr, k*np.ones_like(w, dtype=int)))

    c1 = fits.Column(name='Wavelength', format='D', array=w_arr, unit='Angstroms')
    c2 = fits.Column(name='Flux', format='D', array=f_arr, unit='Counts')
    c3 = fits.Column(name='Order', format='I', array=o_arr)

    tbhdu = fits.BinTableHDU.from_columns([c1,c2,c3])
    tbhdu.writeto(outfile, clobber=True)

def extract_order(ccd, order_frame, n_order, ws, shift_dict, y1=3, y2=10, target=True, interp=False):
    """Given a wavelength solution and offset, extract the order

    """
    hrs = HRSOrder(n_order)
    hrs.set_order_from_array(order_frame.data)
    hrs.set_flux_from_array(ccd.data, flux_unit=ccd.unit)
    hrs.set_target(target)
    data, coef = hrs.create_box(hrs.flux, interp=interp)

    xarr = np.arange(len(data[0]))
    warr = ws(xarr)
    flux = np.zeros_like(xarr, dtype=float)
    weight = 0
    for i in shift_dict.keys():
        if i < len(data) and i >= y1 and i <= y2:
            m = shift_dict[i]
	    shift_flux = np.interp(xarr, m(xarr), data[i])
            data[i] = shift_flux
            flux += shift_flux * np.median(shift_flux)
            weight += np.median(shift_flux)
    pickle.dump(data, open('box_%i.pkl' % n_order, 'w'))
    return warr, flux/weight



#############################33

def extract_order2(ccd, order_frame, n_order, ws, shift_dict, y1=3, y2=10, target=True, interp=False):
    """Given a wavelength solution and offset, extract the order

    """
    hrs = HRSOrder(n_order)
    hrs.set_order_from_array(order_frame.data)
    hrs.set_flux_from_array(ccd.data, flux_unit=ccd.unit)
    hrs.set_target(target)
    data, coef = hrs.create_box(hrs.flux, interp=interp)

    xarr = np.arange(len(data[0]))
    warr = ws(xarr)
    flux = np.zeros_like(xarr, dtype=float)
    weight = 0
    for i in shift_dict.keys():
        if i < len(data) and i >= y1 and i <= y2:
            m = shift_dict[i]
	    shift_flux = np.interp(xarr, m(xarr), data[i])
            data[i] = shift_flux
            flux += shift_flux * np.median(shift_flux)
            weight += np.median(shift_flux)
    pickle.dump(data, open('box_%i.pkl' % n_order, 'w'))
    return xarr, warr, flux/weight



def extract_orders_2D(ccd, order_frame, all_orders, ws, shift_dict, y1=3, y2=10, target=True, interp=False):
    """Given a wavelength solution and offset, extract the order

    """
    xarr_tab=np.array([])
    warr_tab=np.array([])
    order_tab=np.array([])
    flux_tab=np.array([])
    for n_order in all_orders:
        n_order=int(n_order)
        hrs = HRSOrder(n_order)
        hrs.set_order_from_array(order_frame.data)
        hrs.set_flux_from_array(ccd.data, flux_unit=ccd.unit)
        hrs.set_target(target)
        data, coef = hrs.create_box(hrs.flux, interp=interp)

        xarr = np.arange(len(data[0]))
        warr = ws(xarr)
        flux = np.zeros_like(xarr, dtype=float)
        weight = 0
        for i in shift_dict.keys():
            if i < len(data) and i >= y1 and i <= y2:
                m = shift_dict[i]
	        shift_flux = np.interp(xarr, m(xarr), data[i])
                data[i] = shift_flux
                flux += shift_flux * np.median(shift_flux)
                weight += np.median(shift_flux)
        pickle.dump(data, open('box_%i.pkl' % n_order, 'w'))
        xarr_tab=np.append(xarr_tab,xarr)
        oarr=np.ones_like(xarr)*float(n_order)
        order_tab=np.append(order_tab,oarr)
        warr_tab=np.append(warr_tab,warr)
        flux_tab=np.append(flux_tab,flux/weight)
 
    p_init = mod.models.Polynomial2D(degree=4)
    fit_p = mod.fitting.LevMarLSQFitter()
    p = fit_p(p_init,xarr_tab, order_tab, warr_tab)
    fitted=p(xarr_tab, order_tab)
 #   for i in range(len(xarr_tab)):
 #       print xarr_tab[i],order_tab[i],warr_tab[i],fitted[i],warr_tab[i]-fitted[i]
    print "ASDF"
    print all_orders,"CCCC"
    pl.figure()
    pl.subplot(1, 3, 1)
    pl.scatter(xarr_tab, order_tab, c=warr_tab,edgecolor="None")
    pl.xlabel("X")
    pl.ylabel("ORDER")
    pl.colorbar(orientation='horizontal')
    pl.title("Data")
    pl.subplot(1, 3, 2)
    pl.scatter(xarr_tab, order_tab, c=p(xarr_tab, order_tab),edgecolor="None")
    pl.xlabel("X")
    pl.ylabel("ORDER")
    pl.colorbar(orientation='horizontal')
    pl.title("Model")
    pl.subplot(1, 3, 3)
    pl.scatter(xarr_tab, order_tab, c=(warr_tab-p(xarr_tab, order_tab)),edgecolor="None")
    pl.xlabel("X")
    pl.ylabel("ORDER")
    pl.colorbar(orientation='horizontal')
    pl.title("Residual")
    pl.show()
    return True



def extract2(ccd, order_frame, soldir, target='upper', interp=False):
    """Extract all of the orders and create a spectra table

    """
    if target=='upper': 
       target=True
    else:
       target=False

    #set up the orders
    min_order = int(order_frame.data[order_frame.data>0].min())
    max_order = int(order_frame.data[order_frame.data>0].max())
    xarr_tab=np.array([])
    warr_tab=np.array([])
    order_tab=np.array([])
    flux_tab=np.array([])
    sp_dict = {}
    for n_order in np.arange(min_order, max_order):
        try:
            shift_dict, ws = pickle.load(open(soldir+'sol_%i.pkl' % n_order))
        except:
            continue
        x, w, f = extract_order2(ccd, order_frame, n_order, ws, shift_dict, target=target, interp=interp)
        xarr_tab=np.append(xarr_tab,x)
        oarr=np.ones_like(x)*float(n_order)
        order_tab=np.append(order_tab,oarr)
        warr_tab=np.append(warr_tab,w)
        flux_tab=np.append(flux_tab,f)
	sp_dict[n_order] = [w,f]

    p_init = mod.models.Polynomial2D(degree=5)
    fit_p = mod.fitting.LevMarLSQFitter()
#    fit_p = mod.fitting.LinearLSQFitter()
    p = fit_p(p_init,xarr_tab, order_tab, warr_tab)

    residuals=np.abs(warr_tab-p(xarr_tab, order_tab))
    iterations=2
    for tmp in range(iterations):
        std=np.std(residuals)
        mean=np.mean(residuals)
        print std,mean
        i_iter=np.where(residuals < mean+3.*std)
        p = fit_p(p_init,xarr_tab[i_iter], order_tab[i_iter], warr_tab[i_iter])
        residuals=np.abs(warr_tab[i_iter]-p(xarr_tab[i_iter], order_tab[i_iter]))
    
#    for i in range(len(xarr_tab)):
#        print xarr_tab[i],order_tab[i],warr_tab[i],fitted[i],warr_tab[i]-fitted[i]
    pl.figure()
    pl.subplot(1, 3, 1)
    pl.scatter(xarr_tab, order_tab, c=warr_tab,edgecolor="None")
    pl.xlabel("X")
    pl.ylabel("ORDER")
    pl.colorbar(orientation='horizontal')
    pl.title("Single fit")
    pl.subplot(1, 3, 2)
    pl.scatter(xarr_tab, order_tab, c=p(xarr_tab, order_tab),edgecolor="None")
    pl.xlabel("X")
    pl.colorbar(orientation='horizontal')
    pl.title("2D Fit")
    pl.subplot(1, 3, 3)
    pl.scatter(xarr_tab, order_tab, c=(warr_tab-p(xarr_tab, order_tab)),edgecolor="None")
    pl.xlabel("X")
    pl.colorbar(orientation='horizontal')
    pl.title("Residual")
    pl.show()

#    extract_orders_2D(ccd, order_frame, np.arange(min_order, max_order), ws, shift_dict, target=target, interp=interp)
    return sp_dict



####################3

def extract(ccd, order_frame, soldir, target='upper', interp=False):
    """Extract all of the orders and create a spectra table

    """
    if target=='upper': 
       target=True
    else:
       target=False

    #set up the orders
    min_order = int(order_frame.data[order_frame.data>0].min())
    max_order = int(order_frame.data[order_frame.data>0].max())
    sp_dict = {}
    for n_order in np.arange(min_order, max_order):
        try:
            shift_dict, ws = pickle.load(open(soldir+'sol_%i.pkl' % n_order))
        except:
            continue
        w, f = extract_order(ccd, order_frame, n_order, ws, shift_dict, target=target, interp=interp)
	sp_dict[n_order] = [w,f]

    extract_orders_2D(ccd, order_frame, np.arange(min_order, max_order), ws, shift_dict, target=target, interp=interp)
    return sp_dict


if __name__=='__main__':

   
    ccd = CCDData.read(sys.argv[1])
    order_frame = CCDData.read(sys.argv[2], unit=u.adu)
    soldir = sys.argv[3]

    rm, xpos, target, res, w_c, y1, y2 =  mode_setup_information(ccd.header)
    sp_dict = extract2(ccd, order_frame, soldir, interp=True, target=target)
    outfile = sys.argv[1].replace('.fits', '_spec.fits')

    write_spdict(outfile, sp_dict)

