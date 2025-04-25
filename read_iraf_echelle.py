#!/usr/bin/env python3
from sys import argv, exit
from astropy.io import fits
import numpy as np
import warnings
warnings.filterwarnings("ignore")

import matplotlib.pyplot as plt
import matplotlib as mpl

fontsize = 8
mpl.rcParams['xtick.labelsize'] = fontsize
mpl.rcParams['ytick.labelsize'] = fontsize
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'cm'
mpl.rcParams['text.usetex'] = True

def read_multispec(input_file):
    """
    This function reads the input file and returns wavelength and flux.
    Recognize IRAF multipspec spectra with different types of dispersion solution
    ES. 2020-10-21
    """
    try:
        hdu = fits.open(input_file)
    except Exception:
        print("Error while opening the input file")
    finally:
        header = hdu[0].header
        spectrum = hdu[0].data
        hdu.close()
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
    return waves,spectrum

def nonlinearwave(nwave, specstr):
    """
    This function is a modified version of the corresponding unit from readmultispec.py
    (https://raw.githubusercontent.com/kgullikson88/General/master/readmultispec.py)
    Eugene Semenko, 2020-10-21

    Compute non-linear wavelengths from multispec string
    Returns wavelength array and dispersion fields.
    Raises a ValueError if it can't understand the dispersion string.
    """

    fields = specstr.split()
    if int(fields[2]) != 2:
        raise ValueError('Not nonlinear dispersion: dtype=' + fields[2])
    if len(fields) < 12:
        raise ValueError('Bad spectrum format (only %d fields)' % len(fields))
    wt = float(fields[9])
    w0 = float(fields[10])
    ftype = int(fields[11])
    if ftype == 3:
        # cubic spline
        if len(fields) < 15:
            raise ValueError('Bad spline format (only %d fields)' % len(fields))
        npieces = int(fields[12])
        pmin = float(fields[13])
        pmax = float(fields[14])
        if len(fields) != 15 + npieces + 3:
            raise ValueError('Bad order-%d spline format (%d fields)' % (npieces, len(fields)))
        coeff = np.asarray(fields[15:], dtype=float)
        # normalized x coordinates
        s = (np.arange(nwave, dtype=float) + 1 - pmin) / (pmax - pmin) * npieces
        j = s.astype(int).clip(0, npieces - 1)
        a = (j + 1) - s
        b = s - j
        x0 = a ** 3
        x1 = 1 + 3 * a * (1 + a * b)
        x2 = 1 + 3 * b * (1 + a * b)
        x3 = b ** 3
        wave = coeff[j] * x0 + coeff[j + 1] * x1 + coeff[j + 2] * x2 + coeff[j + 3] * x3
    elif ftype == 1 or ftype == 2:
        # chebyshev or legendre polynomial
        # legendre not tested yet
        if len(fields) < 15:
            raise ValueError('Bad polynomial format (only %d fields)' % len(fields))
        order = int(fields[12])
        pmin = float(fields[13])
        pmax = float(fields[14])
        if len(fields) != 15 + order:
            # raise ValueError('Bad order-%d polynomial format (%d fields)' % (order, len(fields)))
            order = len(fields) - 15
        coeff = np.asarray(fields[15:], dtype=float)
        # normalized x coordinates
        pmiddle = (pmax + pmin) / 2
        prange = pmax - pmin
        x = (np.arange(nwave, dtype=float) + 1 - pmiddle) / (prange / 2)
        p0 = np.ones(nwave, dtype=float)
        p1 = x
        wave = p0 * coeff[0] + p1 * coeff[1]
        for i in range(2, order):
            if ftype == 1:
                # chebyshev
                p2 = 2 * x * p1 - p0
            else:
                # legendre
                p2 = ((2 * i - 1) * x * p1 - (i - 1) * p0) / i
            wave = wave + p2 * coeff[i]
            p0 = p1
            p1 = p2
    else:
        raise ValueError('Cannot handle dispersion function of type %d' % ftype)
    return wave

def plot_order(w, r, ordnum):
    nord = np.shape(w)[0]
    fig = plt.figure(figsize=(15,3), tight_layout=True)
    ax = fig.add_subplot(1,1,1)
    ax.set_xlabel(r"Wavelength [\AA]")
    ax.set_ylabel("Intensity")
    if ordnum == -99:
        for i in range(nord):
            ax.plot(w[i], r[i], lw=0.9, ls='-')
    elif ordnum >= 0 and ordnum <= nord-1:
        ax.plot(w[ordnum], r[ordnum], ls='-', lw=0.9, color='red')
    else:
        print("Wrong number of the order")
    plt.show()
    return None

if __name__ == "__main__":
    if len(argv) >= 3:
        wl, sp = read_multispec(argv[1])
        plot_order(wl, sp, int(argv[2]))
    exit(0)
