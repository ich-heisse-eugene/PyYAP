#!/usr/bin/env python3
from sys import argv, exit
import os
from astropy.io import fits
import numpy as np
from numpy.polynomial.legendre import legfit, legval
from numpy.polynomial.chebyshev import chebfit, chebval

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec

import warnings
warnings.filterwarnings("ignore")

fontsize = 7
mpl.rcParams['xtick.labelsize'] = fontsize
mpl.rcParams['ytick.labelsize'] = fontsize

def read_multispec(input_file):
    """
    This function reads the input file and returns wavelength and flux.
    Recognizes IRAF multipspec spectra with different types of dispersion solution
    ES. 2020-10-21
    """
    with fits.open(input_file) as hdu:
        header = hdu[0].header
        spectrum = hdu[0].data

    sizes = np.shape(spectrum)
    if len(sizes) == 1:
        nspec = 1
        npix = sizes[0]
    elif len(sizes) == 2:
        nspec = sizes[0]
        npix = sizes[1]
    elif len(sizes) >=3:
        nspec = sizes[-2]
        npix = sizes[-1]
        spectrum = spectrum[0]
    waves = np.zeros((nspec, npix), dtype=float)
    # try to recognize the type of dispersion
    if 'CTYPE1' in header:
        if header['CTYPE1'].strip() == 'LINEAR': # Linear dispersion solution
            crpix1 = header['CRPIX1']
            crval1 = header['CRVAL1']
            cd1_1 = header['CD1_1']
            wave = (np.arange(npix, dtype=float) + 1 - crpix1) * cd1_1 + crval1
            for i in range(nspec):
                waves[i, :] = wave
            if 'DC-FLAG' in header:
                if header['DC-FLAG'] == 1:
                    waves = 10**waves
        elif header['CTYPE1'].strip() == 'MULTISPE': # IRAF multispec data
            try:
                wat2 = header['WAT2_*']
            except Exception:
                print("Header does not contain keywords required for multispec data")
            finally:
                count_keys = len(wat2)
            long_wat2 = ""
            wave_params = np.zeros((nspec, 24), dtype=float)
            for i in wat2:
                key = header[i].replace('\'', '')
                if len(key) < 68: key += ' '
                long_wat2 += key
            for i in range(nspec):
                idx_b = long_wat2.find("\"", long_wat2.find("spec"+str(i+1)+" ="), -1)
                idx_e = long_wat2.find("\"", idx_b+1, -1)
                temparr = np.asarray(long_wat2[idx_b+1:idx_e].split())
                wave_params[i, 0:len(temparr)] = temparr
                if wave_params[i, 2] == 0 or wave_params[i, 2] == 1:
                    waves[i, :] = np.arange(npix, dtype=float) * wave_params[i, 4] \
                                + wave_params[i, 3]
                    if wave_params[i, 2] == 1:
                        waves[i, :] = 10**waves[i, :]
                else: # Non-linear solution. Not tested
                    waves[i, :] = nonlinearwave(npix, long_wat2[idx_b+1:idx_e])
        elif header['CTYPE1'].strip() == 'PIXEL':
            waves[:,:] = np.arange(npix)+1
    return waves,spectrum,header

def fit_poly(w, r, type, order):
    if type == "legendre":
        coef = legfit(w, r, order)
        return coef
    elif type == "chebyshev":
        coef = chebfit(w, r, order)
        return coef

def fit_cont(w, type, coef):
    if type == "legendre":
        return legval(w, coef)
    elif type == "chebyshev":
        return chebval(w, coef)

def reject_points(w, r, cont, low_rej, high_rej, func, ord):
    resid = r - cont
    stdr = np.std(resid)
    idx = np.where((resid >= -low_rej * stdr) & (resid <= high_rej * stdr))
    coef = fit_poly(w[idx], r[idx], func, ord)
    return w[idx], r[idx], cont, coef

def normalise_flat(source_file, fit_func, fit_ord, fit_niter, fit_low_rej, fit_high_rej):
    w_init, r_init, hdr_init = read_multispec(source_file)
    cont_lev = np.zeros(np.shape(w_init), dtype=r_init[0].dtype)
    nord = np.shape(w_init)[0]
    for ord in range(nord):
        w_tmp = w_init[ord, :]
        r_tmp = r_init[ord, :]
        for j in range(fit_niter+1):
            coef = fit_poly(w_tmp, r_tmp, fit_func, fit_ord)
            cont = fit_cont(w_tmp, fit_func, coef)
            w_tmp, r_tmp, cont, coef = reject_points(w_tmp, r_tmp, cont, fit_low_rej, fit_high_rej, fit_func, fit_ord)
            cont_cur = fit_cont(w_tmp, fit_func, coef)
            idx_wrong = np.where(cont_cur == 0.)
            cont_cur[idx_wrong] = r_tmp[idx_wrong]
        cont_lev[ord, :] = fit_cont(w_init[ord, :], fit_func, coef)
        # plt.plot(w_init[ord, :], r_init[ord, :]/cont_lev[ord, :], 'b-', lw=0.7)
        # plt.plot(w_init[ord, :], cont_lev[ord, :], 'r-', lw=1.7)
        # plt.show()
    return cont_lev/r_init

def remove_blz(file_spec, file_blaze, file_corr):
    _, r_s, hdr_s = read_multispec(file_spec)
    _, r_b, _ = read_multispec(file_blaze)

    r_cor = r_s / r_b
    hdr_s['HISTORY'] = 'Blaze normalization of the spectrum'
    hdu = fits.PrimaryHDU(r_cor)
    hdu.header = hdr_s.copy()
    hdu.writeto(file_corr, overwrite=True)
    return file_corr

if __name__ == "__main__":
    # norm = normalise_flat(argv[1], "chebyshev", 12, 0, 3., 3.)
    wf, rf, hdrf = read_multispec(argv[1])
    w, r, hdr = read_multispec(argv[2])
    nord = np.shape(w)[0]
    for ord in range(nord):
        plt.plot(w[ord, :], r[ord, :]/rf[ord, :], lw=0.5)
    plt.show()
    # hdu = fits.PrimaryHDU(r/rf)
    # hdu.header = hdr.copy()
    # hdu.writeto("flat_norm_"+argv[2], overwrite=True)
    exit(0)

