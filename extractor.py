import astropy.io.fits as pyfits
import os
import numpy as np
import scipy.interpolate
import shutil

import logging

import matplotlib
import matplotlib.pyplot as plt

from skimage.exposure import equalize_hist
from numpy.polynomial.chebyshev import chebval

import warnings
warnings.simplefilter("ignore")

##################################################################
def fox(file_name, sflat_image, ap_file, ex_type, aperture):

    print(f"Extraction spectra from {file_name} started")
    logging.info(f"Extraction spectra from {file_name} started")

    #read file with super flat
    hdulist = pyfits.open(sflat_image)
    s_flat_data = hdulist[0].data.copy()
    prihdr = hdulist[0].header
    if 'RDNOISE' in prihdr:
        RN = prihdr['RDNOISE']
        print("Read noise = ", RN)
        logging.info(f"Read noise = {RN:.1f}")
    else:
        print("Read noise not found in", sflat_image, " header.")
        logging.warning(f"Read noise not found in {sflat_image} header.")
        return (None)
    if 'GAIN' in prihdr:
        gain = prihdr['GAIN']
        print("Gain = ", gain)
        logging.info(f"Gain = {gain:.1f}")
    else:
        print("Gain not found in", sflat_image, " header.")
        logging.warning(f"Gain not found in {sflat_image} header.")
        return (None)
    hdulist.close()

    #read file with spectra
    hdulist = pyfits.open(file_name)
    spectra_data = hdulist[0].data.copy()
    prihdr = hdulist[0].header
    hdulist.close()

    if s_flat_data.shape[0]!=spectra_data.shape[0] or s_flat_data.shape[1]!=spectra_data.shape[1]:
        return (None)

    orders = []
    width = []
    X = []
    Y = []
    #read apertures file
    ii=0
    with open(ap_file, 'r') as f:
        for line in f:
            one_str = line.rsplit('\t')
            if one_str[0]=='Order':
                if ii!=0:
                    X=np.array(X)
                    Y=np.array(Y)
                    orders.append(X)
                    width.append(Y)
                    X = []
                    Y = []
                ii=ii+1
            if one_str[0]!='Order':
                X.append(float(one_str[0]))
                Y.append(float(one_str[1]))
    X=np.array(X)
    Y=np.array(Y)
    orders.append(X)
    width.append(Y)
    X = []
    Y = []
    f.close()
    print(ii, "orders read from file")
    orders = np.array(orders) #array with apretures
    width = np.array(width)

    s_flat_data = np.clip(s_flat_data, 0, 65000)
    spectra_data = np.clip(spectra_data, 0, 65000)

    #flat-relative optimal extraction
    if ex_type=='FOX':
        sq_RN = np.zeros_like(s_flat_data)
        sq_RN.fill((RN/gain)**2)##################################
        sq_RN = sq_RN + spectra_data    #
        weight = np.ones_like(s_flat_data)
        weight = np.divide(weight, sq_RN)
        #multiply arays
        up_sum = np.multiply(weight,s_flat_data)
        up_sum = np.multiply(up_sum,spectra_data)
        down_sum = np.multiply(weight,s_flat_data)
        down_sum = np.multiply(down_sum,s_flat_data)

        #FOX extract orders
        s_ex = np.zeros((orders.shape[0], spectra_data.shape[1]))
        err_ex = np.zeros((orders.shape[0], spectra_data.shape[1]))
        for ii in range(0,s_ex.shape[1]):
            for jj in range(0,s_ex.shape[0]):
                FWHM = width[jj,ii]
                center = orders[jj,ii]
                up_f = center-aperture*FWHM
                down_f = center+aperture*FWHM
                up=round(up_f)+1
                down=round(down_f)-1
                up_part = (up - up_f) - 0.5
                down_part = (down_f - down) - 0.5
                try:
                    _up =   np.sum(  up_sum[int(up):int(down+1),ii])+  up_sum[int(up-1),ii]*up_part+  up_sum[int(down+1),ii]*down_part
                    _down = np.sum(down_sum[int(up):int(down+1),ii])+down_sum[int(up-1),ii]*up_part+down_sum[int(down+1),ii]*down_part
                    s_ex[jj,ii] = _up / _down
                    err_ex[jj,ii] = _up / np.sqrt(_down)
                except:
                    s_ex[jj,ii] = 0
                    err_ex[jj,ii] = 0

    #simple aperture extraction
    if ex_type=='APEX':
        s_ex = np.zeros((orders.shape[0], spectra_data.shape[1]), float)
        err_ex = np.zeros((orders.shape[0], spectra_data.shape[1]))
        for ii in range(0,s_ex.shape[1]):
            for jj in range(0,s_ex.shape[0]):
                FWHM = width[jj,ii]
                center = orders[jj,ii]
                up_f = center-aperture*FWHM
                down_f = center+aperture*FWHM
                up=round(up_f)+1
                down=round(down_f)-1
                up_part = (up - up_f) - 0.5
                down_part = (down_f - down) - 0.5
                try:
                    s_ex[jj,ii] = np.sum(spectra_data[int(up):int(down+1),ii])+spectra_data[int(up-1),ii]*up_part+spectra_data[int(down+1),ii]*down_part
                    err_ex[jj,ii] = s_ex[jj,ii]/np.sqrt(s_ex[jj,ii]+(down_f-up_f)*(RN/gain)**2)
                except Exception as e:
                    print(f"File: {file_name}, Exception: {e}")
                    s_ex[jj,ii] = 0
                    err_ex[jj,ii] = 0

    #PSF-weighted extraction
    if ex_type=='PSFEX':
        s_ex = np.zeros((orders.shape[0], spectra_data.shape[1]), float)
        err_ex = np.zeros((orders.shape[0], spectra_data.shape[1]))
        for ii in range(0,s_ex.shape[1]):#1000
            for jj in range(0,s_ex.shape[0]):#29
                FWHM = width[jj,ii]
                center = orders[jj,ii]
                sigma = FWHM/2.35482    #gauss model constant

                up_f   = center - aperture*FWHM
                down_f = center + aperture*FWHM
                up=round(up_f)
                down=round(down_f)
                up_part = 0.5+(up-up_f)
                down_part = 0.5+(down_f-down)
                try:
                    ROI = spectra_data[int(up):int(down+1),ii]
                    X = np.arange(0,len(ROI))
                    x0=center-up
                    profile = np.exp(-1*(X-x0)*(X-x0)/(2*sigma*sigma))     #gauss model
                    norm_sum = profile[0]*up_part + profile[len(ROI)-1]*down_part + np.sum(profile[1:len(ROI)-1])
                    norm_profile =profile / norm_sum

                    variance = np.zeros_like(ROI)
                    variance.fill((RN/gain)**2)
                    variance = variance + ROI
                    down_arr = np.divide(np.multiply(norm_profile, norm_profile), variance)
                    up_arr   = np.divide(np.multiply(ROI,          norm_profile), variance)

                    ##extraction
                    up_sum =     up_arr[0]*up_part +   up_arr[len(ROI)-1]*down_part + np.sum(  up_arr[1:len(ROI)-1])
                    down_sum = down_arr[0]*up_part + down_arr[len(ROI)-1]*down_part + np.sum(down_arr[1:len(ROI)-1])
                    s_ex[jj,ii] = up_sum / down_sum
                    ##error
                    err_arr = norm_profile*norm_profile/variance
                    err_sum = err_arr[0]*up_part + err_arr[len(ROI)-1]*down_part + np.sum(err_arr[1:len(ROI)-1])
                    err_arr = (norm_profile/variance)/err_sum
                    err_arr = err_arr*err_arr*variance
                    err_sum = err_arr[0]*up_part + err_arr[len(ROI)-1]*down_part + np.sum(err_arr[1:len(ROI)-1])
                    err_ex[jj,ii] = s_ex[jj,ii] / np.sqrt(err_sum)
                except:
                    s_ex[jj,ii] = 0
                    err_ex[jj,ii] = 0

    s_ex=np.float32(s_ex)
    s_ex = np.fliplr(s_ex)


    hdu = pyfits.PrimaryHDU(s_ex)
    hdu.header = prihdr
    hdu.header['HISTORY'] = 'extracted by '+ex_type
    hdu.header['HISTORY'] = 'flip X'
    for ii in range (0, s_ex.shape[0]):
        hdu.header['APNUM'+str(ii+1)] = str(ii+1)+' '+str(ii+1)+' '+str(np.min(orders[ii,:]))+' '+str(np.max(orders[ii,:]))

    hdu.header['WCSDIM'] = 2
    hdu.header['CTYPE1'] = 'PIXEL'
    hdu.header['CTYPE2'] = 'LINEAR'
    hdu.header['CRPIX1'] = 1.0
    hdu.header['CDELT1'] = 1.0
    hdu.header['CDELT2'] = 1.0
    hdu.header['CD1_1'] = 1.0
    hdu.header['CD2_2'] = 1.0
    hdu.header['LTM1_1'] = 1.0
    hdu.header['LTM2_2'] = 1.0
    hdu.header['WAT0_001'] = 'system=equispec'
    hdu.header['WAT1_001'] = 'wtype=linear label=Pixel'
    hdu.header['WAT2_001'] = 'wtype=linear'

    new_file = os.path.splitext(file_name)[0]+ '_ec.fits'
    hdu.writeto(new_file, overwrite=True)

    err_ex=np.float32(err_ex)
    err_ex = np.fliplr(err_ex)
    hdu = pyfits.PrimaryHDU(err_ex)
    hdu.header = prihdr
    hdu.header['IMAGETYP'] = 'S/N map'
    err_file = os.path.splitext(file_name)[0] + '_err.fits'
    hdu.writeto(err_file, overwrite=True)

    return(new_file, err_file)
