from astropy.io import fits
import os
import numpy as np
import scipy.interpolate
import shutil
import multiprocessing as mp
import logging

import matplotlib
import matplotlib.pyplot as plt

from skimage.exposure import equalize_hist

import warnings
warnings.simplefilter("ignore")

##################################################################

def read_traces(x_coo, ap_file):
    Y = []
    FWHM = []
    with open(ap_file, 'r') as fp:
        f = fp.readlines()
        p = f[3].strip().rsplit()
        if len(p) == 4:
            n_orders = int(p[0]); poly_order = int(p[1])
            adaptive = p[2]; ap_size = float(p[3])
        else:
            print("Problem of reading the file with traces")
        if adaptive == 'True':
            for i in range(4, 4+n_orders):
                p = f[i].strip().rsplit()
                poly_trace_coef = np.asarray(p[3:3+poly_order+1], dtype=float)
                poly_width_coef = np.asarray(p[3+poly_order+1:], dtype=float)
                Y.append(np.polyval(poly_trace_coef, x_coo))
                FWHM.append(np.polyval(poly_width_coef, x_coo))
        else:
            for i in range(4, 4+n_orders):
                p = f[i].strip().rsplit()
                poly_trace_coef = np.asarray(p[3:3+poly_order+1], dtype=float)
                Y.append(np.polyval(poly_trace_coef, x_coo))
                FWHM.append(np.repeat(float(p[3+poly_order+1]), len(x_coo)))
    return np.asarray(Y), np.asarray(FWHM)

def fox(file_name, sflat_image, ap_file, ex_type, aperture, queue):
    with fits.open(sflat_image) as hdulist:
        s_flat_data = hdulist[0].data.copy()
        prihdr = hdulist[0].header
        if 'RDNOISE' in prihdr:
            RN = prihdr['RDNOISE']
        else:
            RN = -999
            return None
        if 'GAIN' in prihdr:
            gain = prihdr['GAIN']
        else:
            gain = -999
            return None

    #read file with spectra
    with fits.open(file_name) as hdulist:
        spectra_data = hdulist[0].data.copy()
        prihdr = hdulist[0].header

    if s_flat_data.shape[0]!=spectra_data.shape[0] or s_flat_data.shape[1]!=spectra_data.shape[1]:
        return (None)

    orders, width = read_traces(np.arange(s_flat_data.shape[1]), ap_file)
    print(f"Extraction of {len(orders)} spectral orders from {file_name} started. Gain = {gain:.1f}\tReadout noise = {RN:.2f}")
    queue.put((logging.INFO, f"Extraction of {len(orders)} spectral orders from {file_name} started. Gain = {gain:.1f}\tReadout noise = {RN:.2f}"))

    s_flat_data = np.clip(s_flat_data, 0, 65000)
    spectra_data = np.clip(spectra_data, 0, 65000)

    #flat-relative optimal extraction
    if ex_type=='FOX':
        sq_RN = np.zeros_like(s_flat_data)
        sq_RN.fill((RN/gain)**2)
        sq_RN = sq_RN + spectra_data
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

    hdu = fits.PrimaryHDU(s_ex)
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
    hdu = fits.PrimaryHDU(err_ex)
    hdu.header = prihdr
    hdu.header['IMAGETYP'] = 'S/N map'
    err_file = os.path.splitext(file_name)[0] + '_err.fits'
    hdu.writeto(err_file, overwrite=True)
    return (new_file, err_file)

def extract_multi(Path2Data, Path2Temp, list_name, out_list_name, out_err_list_name, conf, flat_name, queue):
    ap_file = os.path.join(Path2Temp, 'traces.txt')
    aperture = float(conf['aperture'])
    ex_type = conf['ex_type']
    view = eval(conf['view'])
    out_list = []
    with open(list_name, 'r') as f:
        proc_args = [(Path2Temp, line.strip(), flat_name, ap_file, ex_type, aperture, queue) for line in f]
        nCPUs = os.cpu_count()
        if 'threading' in conf and eval(conf['threading']) and nCPUs > 2:
            with mp.Pool(processes=nCPUs) as pool:
                res_async = pool.starmap_async(process_multi, proc_args, chunksize=nCPUs)
                res_async.wait()
                out_list.extend(res_async.get())
        else:
            for item in proc_args:
                res_mono = process_multi(Path2Temp, item[1], flat_name, ap_file, ex_type, aperture, queue)
                out_list.extend([res_mono])
    f_out = open(out_list_name, 'a')
    fe_out = open(out_err_list_name, 'a')
    for item in out_list:
        print(item[0], file=f_out)
        print(item[1], file=fe_out)
    f_out.close()
    fe_out.close()
    os.remove(list_name)
    print("Extraction complete")
    queue.put((logging.INFO, "Extraction complete"))
    return None

def process_multi(Path2Temp, name, flat_name, ap_file, ex_type, aperture, queue):
    status, status_err = fox(name, flat_name, ap_file, ex_type, aperture, queue)
    print(f"Extracted spectrum saved in {status}")
    queue.put((logging.INFO, f"Extracted spectrum saved in {status}"))
    shutil.move(os.fspath(name), os.fspath(Path2Temp))
    return [status, status_err]
