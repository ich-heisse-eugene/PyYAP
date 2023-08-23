"""
    Programme estimates the resolving power from fitting of ThAr lines and produces
    a PDF file with report
    Written by Eugene Semenko
    Last modification: 2021-01-22
    Changelog:
    2020-12-15. Fixed normalisation of an average profile in each order. Added
                labels to the axis. ES
    2021-01-22. Replaced some symbols with letters in Latin script.
"""
import datetime
from astropy.io import fits as pyfits
import numpy as np
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import warnings
warnings.filterwarnings("ignore")

import logging

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages

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
        hdu = pyfits.open(input_file)
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

def resol_in_order(wave, spec, n):
    """
        Computes resolution in a separate order
    """
    half_win = 5
    pixarr = np.arange(len(wave), dtype=float)
    idx = np.where(spec < np.mean(spec))
    if len(idx[0]) == 0:
        return -1,-1,-1,-1,-1
    try:
        coef = np.polyfit(wave[idx], spec[idx], 4) # bkg
    except:
        return -1,-1,-1,-1,-1
    else:
        signal = spec - np.polyval(coef, wave) # with subtracted bkg
    lind,_ = find_peaks(signal/np.max(signal), height=0.05, distance=10)
    res = []
    wl = []
    dwl = []
    x = []
    y = []
    for i in lind:
        if i - half_win > 0 and i + half_win < pixarr[-1]:
            res_cur, wl0, dwl_cur, x_cur, y_cur = fit_line(wave[i-half_win : i+half_win+1], \
                                                           spec[i-half_win : i+half_win+1])
            if res_cur != -1:
                res.append(res_cur)
                wl.append(wl0)
                dwl.append(dwl_cur)
                x.append(x_cur)
                y.append(y_cur)
    if len(res) != 0:
        print(f"Order #{n:2.0f} {wave[0]:.1f}-{wave[-1]:.1f}A: Median R = {np.median(res):.0f}, dl = {np.median(dwl):.4f}. {len(res)} lines measured")
        logging.info(f"Order #{n:2.0f} {wave[0]:.1f}-{wave[-1]:.1f}A: Median R = {np.median(res):.0f}, dl = {np.median(dwl):.4f}. {len(res)} lines measured")
    else:
        return -1,-1,-1,-1,-1
    return res, wl, dwl, x, y

def fit_line(w, y):
    gauss = lambda x,a,b,c,d: a*np.exp(-(x - b)**2/(2. * c**2)) + d
    dlam = (w[-1] - w[0]) / len(w)
    x_out = []
    y_out = []
    try:
        p0 = np.array([np.max(y), np.mean(w), 3.*dlam, 1.])
        popt, pcov = curve_fit(gauss, w, y, p0, maxfev=10000)
    except Exception:
        return -1, -1, -1, x_out, y_out
    else:
        fwhm = 2*np.sqrt(2*np.log(2))*popt[2]
        if fwhm >= 2*dlam and fwhm <= 6*dlam:
            res_fwhm = int(popt[1] / fwhm)
            x_cen = (w - popt[1]) / dlam
            y = (y - popt[3]) / popt[0]
            x_out.append(x_cen)
            y_out.append(y)
        else:
            return -1, -1, -1, x_out, y_out
    return res_fwhm, popt[0], dlam, x_out, y_out

def make_report_resol(w, sp, input_file, view=False):
    print("Test of spectral resolution:")
    logging.info("Test of spectral resolution:")
    orders = np.arange(np.shape(w)[0], dtype=int)
    norders = len(orders)
    res = np.zeros(norders)
    dwl = np.zeros(norders)
    gauss = lambda x,a,b,c,d: a*np.exp(-((x - b)/(2.*c))**2) + d
    result_wl0 = np.array([])
    result_dwl0 = np.array([])
    result_x = np.array([])
    result_y = np.array([])
    result_resol = np.array([])
    nlines = 0
    # Multipage PDF output
    with PdfPages(input_file.replace(".fits", ".report.pdf")) as pdf:
        d = pdf.infodict()
        d['Title'] = f"Spectral resolution measured from file {input_file}"
        d['ModDate'] = datetime.datetime.today()
        nx = 2 # Number of columns in output
        ny = 5 # number of rows in output
        nplots = nx*ny # number of plots per page
        npages = int(np.ceil(norders / nplots))
        cur_order = 0
        for page in range(npages):
            fig = plt.figure(constrained_layout=True, dpi=300, tight_layout=True)
            fig.set_size_inches(8.27, 11.69, forward=True)
            spec = gridspec.GridSpec(ncols=nx, nrows=ny, figure=fig, hspace=0.1)
            cur_row = 0; cur_col = 0
            for i in orders[page+(page*nplots): page+(page*nplots)+nplots]:
                res_ord, wl_ord, dwl_ord, x_ord, y_ord = resol_in_order(w[i,:], sp[i,:], i)
                if res_ord != -1:
                    x_ord = np.array(x_ord).flatten()
                    y_ord = np.array(y_ord).flatten()
                    result_resol = np.hstack((result_resol, res_ord))
                    result_x = np.hstack((result_x, x_ord))
                    result_y = np.hstack((result_y, y_ord))
                    result_wl0 = np.hstack((result_wl0, np.array(wl_ord).flatten()))
                    result_dwl0 = np.hstack((result_dwl0, np.array(dwl_ord).flatten()))
                    nlines += len(res_ord)
                    f_ax = fig.add_subplot(spec[cur_row, cur_col])
                    # print(f"Order {i}, R={np.mean(res_ord):.0f}")
                    if len(x_ord) != 0 and len(y_ord) != 0:
                        f_ax.plot(x_ord, y_ord, 'b.', ms=0.7)
                    f_ax.set_title(f"Order {i}, {w[i,0]:.0f}-{w[i,-1]:.0f}Å", fontsize=fontsize)
                    f_ax.set_xlabel("Pixel", fontsize=fontsize)
                    f_ax.set_ylabel("Normalized intensity", fontsize=fontsize)
                    try:
                        p0 = np.array([1, 0, 1., 1.])
                        popt, pcov = curve_fit(gauss, x_ord, y_ord, p0, maxfev=10000)
                    except RuntimeError or ValueError or OptimizeWarning:
                        pass
                    else:
                        fwhm = 2*np.sqrt(2*np.log2(2))*popt[2]
                        xx = np.linspace(np.min(x_ord), np.max(x_ord), 100)
                        f_ax.plot(xx, gauss(xx, *popt), 'r-', lw=0.5)
                        R = np.mean(w[i,:]) / (np.median(dwl_ord)*fwhm)
                        label = f"R={R:.0f}\n $\Delta\lambda$={np.median(dwl_ord):.4f} Å/px\nFWHM={fwhm:.2f} px\n{len(dwl_ord)} lines"
                        f_ax.annotate(label, (2.5, 0.7), va='center', fontsize=fontsize-1)
                    cur_col += 1
                    if cur_col > nx-1:
                        cur_col = 0
                        cur_row += 1
            pdf.savefig()
            plt.close()
        try:
            p0 = np.array([1, 0, 1., 1.])
            popt, pcov = curve_fit(gauss, result_x, result_y, p0, maxfev=10000)
        except RuntimeError or ValueError or OptimizeWarning:
            pass
        finally:
            fwhm = 2*np.sqrt(2*np.log2(2))*popt[2]
            xx = np.linspace(np.min(result_x)-1, np.max(result_x)+1, 100)
            fig = plt.figure(constrained_layout=True, dpi=300, tight_layout=True)
            fig.set_size_inches(8.27, 11.69, forward=True)
            ax = fig.add_subplot(1,1,1)
            ax.plot(result_x, result_y, 'b.', ms=0.5)
            ax.plot(xx, gauss(xx, *popt), 'r-', lw=0.5)
            R = np.mean(w) / (np.median(result_dwl0)*fwhm)
            label = f"Median R={np.median(result_resol):.0f}\nMedian $\Delta\lambda$={np.median(result_dwl0):.4f} \
                               Å/px\n R from fit = {R:.0f} \n FWHM={fwhm:.2f} px\n{nlines} lines"
            ax.annotate(label, (2.5, 0.7), va='center', fontsize=fontsize-1)
            ax.set_xlabel("Pixel")
            ax.set_ylabel("Normalized intensity")
        pdf.savefig()
        plt.close()
    return np.median(result_resol)

def get_sp_resolv(input_file):
    w, sp = read_multispec(input_file)
    R = make_report_resol(w, sp, input_file, False)
    try:
        hdulist = pyfits.open(input_file, mode = 'update')
        prihdr = hdulist[0].header
    except IOError:
        print ("Input/output error. File:", name)
    finally:
        prihdr.set('R', R, 'Median spectral resolution of data')
        prihdr['HISTORY'] = 'Spectral resolution was measured from ThAr spectrum'
        hdulist[0].header = prihdr
        hdulist.close()
    return None
