from astropy.io import fits
import os
import numpy as np
import scipy.ndimage
from scipy.optimize import curve_fit, OptimizeWarning
from scipy import optimize
import matplotlib
import matplotlib.pyplot as plt
import logging

import warnings
warnings.simplefilter("ignore", category=OptimizeWarning)

c_const = 299792.458

def write_disp(name, data, OS, X_Order, Y_Order, prihdr):
    f_name = os.path.splitext(name)[0] + "_disp.txt"
    text_file = open(f_name, "w")
    text_file.write(str(OS)+'\n')
    text_file.write(str(X_Order)+'\n')
    text_file.write(str(Y_Order)+'\n')
    for ii in range(0, len(data)):
        text_file.write(str(data[ii])+'\n')
    text_file.close()
    prihdr['WLSLTN'] = (str(f_name.split(os.sep)[-1]), 'wavelength solution')
    return(prihdr)

####################################################################
def write_features(name, features, prihdr):
    f_name = os.path.splitext(name)[0] + "_features.txt"
    text_file = open(f_name, "w")
    for ii in range(0, len(features)):
        print(f"{int(features[ii][0])}\t{features[ii][1]:.3f}\t{features[ii][2]:.4f}\t{features[ii][3]:.3f}", file=text_file)
    text_file.close()
    prihdr['FEATURES'] = (str(f_name.split(os.sep)[-1]), 'ThAr features for WLSLTN')
    return(prihdr)

####################################################################
def write_new_short(spectrum, features, prihdr, dir_name):
    text_file = open(os.path.join(dir_name, 'thar_new.dat'), "w")
    for ii in range(0, len(features)):
        print(str(int(features[ii][0]))+'\t'+str(round(features[ii][1],3))+'\t'+str(round(features[ii][2],4))+'\t'+str(round(features[ii][3],3)), file=text_file)
    text_file.close()
    hdu = fits.PrimaryHDU(spectrum)
    hdu.header = prihdr
    hdu.writeto(os.path.join(dir_name, 's_thar_new.fits'), overwrite=True)

####################################################################
### checker
def WL_checker(WL, o, OS, features):
    for ii in range (0, len(features)):
            if WL ==features[ii][2] and o+OS==features[ii][0]:
                return (1)
    return (0)

####################################################################
###1D fitting/center search
def center(s_order, x_coo, y_coo, bckgrnd):
    fit_area = FWHM
    x_range_min = x_coo-fit_area
    x_range_max = x_coo+fit_area+1
    if x_range_min < 0:
        x_range_min = 0
    if x_range_max > s_order.shape[0]-1:
        x_range_max = s_order.shape[0]-1
    x_range = np.linspace(x_range_min, x_range_max, (x_range_max - x_range_min), dtype=int)
    ROI = s_order[x_range]
    X = np.arange(ROI.shape[0])
    x0 = ROI.shape[0]//2
    #moffat fitting
    #A - amplitude, B,C - coeff of function (width ...), D - background
    # moffat = lambda x, A, B, C, D, x0: A*(1 + ((x-x0)/B)**2)**(-C)+D
    # p0 = np.array([y_coo, 3, 3, bckgrnd, x0])
    # Gaussian fitting
    gauss = lambda x,a,b,c,d: a*np.exp(-(x - b)**2/(2. * c**2)) + d
    p0 = np.array([y_coo, x0, 1, bckgrnd])
    if y_coo >= threshold*bckgrnd:
        try:
            popt, pcov = curve_fit(gauss, X, ROI, p0, maxfev=10000, method='lm')
        except Exception:
            return 0
        # print(p0)
        # plt.plot(X, ROI, ls='', marker='.', ms=4)
        # xx = np.linspace(X[0], X[-1], 500)
        # plt.plot(xx, moffat(xx, *popt), ls='-', lw=0.5)
        # plt.show()
        return (popt[1]-fit_area)
    else:
        return 0

####################################################################
###search max element in 1d array

def get_max(s_order, x_coo):
    x_range = np.linspace(int(x_coo)-FWHM, int(x_coo)+FWHM+1, 2*FWHM+1, dtype=int)
    ROI =  s_order[x_range]
    x_max = int(x_coo)-FWHM+np.argmax(ROI)
    Max = s_order[x_max]
    Min = np.min(ROI)
    return (x_max, Max, Min)

##################################################################
### 2D dispersion function fitting
def solution(C, Y_Order): #2d chebyshev polynom calculation
    shift = C[0]
    C = np.delete(C,0)
    coeff = np.reshape(C,(-1,Y_Order))#7-Yorder
    return lambda X,O:(np.polynomial.chebyshev.chebval2d(X,O,coeff)-shift)/O

def dispers_p(features, X_Order, Y_Order, view):#, fit_params):
    good=True
    while (good):
        local=np.asarray(features)
        O = np.copy(local[:, 0])
        X = np.copy(local[:, 1])
        W = np.copy(local[:, 2])
        params = np.zeros([1+X_Order*Y_Order])
        errorfunction = lambda C: np.ravel(solution(C, Y_Order)(X,O) - W)
        C, success = optimize.leastsq(errorfunction, params, maxfev=10000)
        residual = solution(C, Y_Order)(X,O)
        residual = W-residual
        good=False
        if np.max(np.absolute(residual))>tolerance/2:
            good=True
            index = np.where(np.absolute(residual)>tolerance/2)[0]
            local = np.delete(local, index.tolist(), 0)
            features = local.tolist()

    rms = round(np.std(residual),5)
    points = O.shape[0]
    print(f"Lines used: {points}")
    print(f"RMS(Angstrom): {rms:.4f}")
    print(f"RMS(km/s): {round(np.std((residual/W)*c_const),5):.2f}")
    logging.info(f"Lines used: {points}")
    logging.info(f"RMS(Angstrom): {rms:.4f}")
    logging.info(f"RMS(km/s): {round(np.std((residual/W)*c_const),5):.2f}")

    if view:
        fig = plt.figure(1, figsize=(16, 5))
        ax = fig.add_subplot(111)
        points = O.shape[0]
        rms = round(np.std(residual),5)
        label_data = f"Chebyshev, Xorder={X_Order}, Yorder={Y_Order}, {points} points, RMS={rms} Angstrom ({np.std((residual/W)*c_const):.2f} km/s)"
        plt.cla()
        ax.set_xlim([np.min(X), np.max(X)])
        ax.plot(X, residual, 'go')
        plt.title(label_data)
        plt.ylabel('Angstrom', fontsize=15)
        plt.xlabel('Pixel', fontsize=15)
        plt.show()

    return (C, features, points, rms)

####################################################################
### auto search lines
def add_lines(spectrum, OS, new_features, thar, disp_params, Y_Order):
    for ii in range(0, spectrum.shape[0]):
        old_len = len(new_features)
        row = spectrum[ii,:]
        # smoothing data
        row_S = scipy.ndimage.gaussian_filter(row, int(FWHM/3))
        row_S = row_S / np.max(row_S)
        row_f = []
        for jj in range (int(FWHM*2),row_S.shape[0]-int(FWHM*2)):
            #search for local extremum of smoothed spectra
            if (row_S[jj]-row_S[jj-1]>0 and row_S[jj+1]-row_S[jj]<0):
                # search for a maximum near clicked position
                offset = get_max(row, jj)
                x_coo = offset[0]
                y_coo = offset[1]
                bckgrnd = offset[2]
                if y_coo >= bckgrnd*threshold:
                    PWL = solution(disp_params, Y_Order)(x_coo,ii+OS)
                    loc_thar = (thar - PWL)**2
                    nearest = np.argmin(loc_thar)
                    if np.sqrt(loc_thar[nearest]) < tolerance*2:
                        offset = center(row, x_coo, y_coo, bckgrnd)
                        if offset!=0:
                            x_coo=round(x_coo+offset,3)
                            PWL = solution(disp_params, Y_Order)(x_coo,ii+OS)
                            loc_thar = (thar - PWL)**2
                            nearest = np.argmin(loc_thar)
                            if np.sqrt(loc_thar[nearest])<tolerance:
                                WL = thar[nearest]
                                if WL_checker(WL, ii, OS, new_features)==0:
                                    new_features.append([ii+OS, x_coo, WL, y_coo])
        print(f"{len(new_features) - old_len:.0f} features found in order {ii+OS:.0f}")
        logging.info(f"{len(new_features) - old_len:.0f} features found in order {ii+OS:.0f}")

    return (new_features)

####################################################################
###search shift of order
def search_shift(spectrum, zero):
    shift = []
    half_max_shift = 15
    for i in range(spectrum.shape[0]):
        ccf = np.correlate(spectrum[i], zero[i], mode='full')
        ccf = ccf / ccf.max()
        x = np.arange(len(ccf))
        y_ccf = ccf[len(ccf)//2 - half_max_shift: len(ccf)//2 + half_max_shift+1]
        x_ccf = np.arange(-half_max_shift, half_max_shift+1, 1)
        # plt.plot(x_ccf, y_ccf, marker='.', ms=5, ls='')
        s_shift = np.argmax(y_ccf) - len(y_ccf)//2
        gauss = lambda x,a,b,c,d: a*np.exp(-(x - b)**2/(2. * c**2)) + d
        p0 = np.array([ccf[s_shift], s_shift, 1., 0.])
        try:
            popt, pcov = curve_fit(gauss, x_ccf, y_ccf, p0, maxfev=10000)
            # xx = np.linspace(x_ccf[0], x_ccf[-1], 500)
            # plt.plot(xx, gauss(xx, *popt), ls='-', lw=0.5)
        except Exception:
            pass
        else:
            if popt[1] >= x_ccf[0] and popt[1] <= x_ccf[-1]:
                shift.append(popt[1])
    # plt.show()
    return np.mean(shift)

####################################################################
###reidentify features from short list
def reidentify_features(s_order, OS, zero_features, line, shift, new_features):
    for ii in range(0,len(zero_features)):
        if zero_features[ii][0] == line+OS:
            WL = zero_features[ii][2]
            x_coo = zero_features[ii][1]
            x_coo_new = x_coo+shift
            try:
                offset = get_max(s_order, np.round(x_coo_new))     #search max pixel near click
                x_coo_new = offset[0]
                y_coo_new= offset[1]
                bckgrnd = offset[2]
                offset = center(s_order, x_coo_new, y_coo_new, bckgrnd)
                if offset!=0:
                    x_coo_new = round(x_coo_new+offset,4)
                    new_features.append([line+OS, x_coo_new, WL, y_coo_new])
            except Exception:
                pass
    return(new_features)

######################################################################
def first_ident(spectrum, zero, zero_features, OS):
    shift = search_shift(spectrum, zero)  #get shift
    print(f"Average shift betweet ThAr epochs is {shift:.2f} pix")
    logging.info(f"Average shift betweet ThAr epochs is {shift:.2f} pix")
    new_features=[]
    #for each order
    found=0
    for ii in range(0, spectrum.shape[0]):
        s_order = spectrum[ii]
        s_order = s_order / np.max(s_order)
        z_order = zero[ii]
        z_order = z_order / np.max(z_order)
        new_features = reidentify_features(s_order, OS, zero_features, ii, shift, new_features)
        print(f"{len(new_features) - found:.0f} features identified in order {ii+OS:.0f}")
        logging.info(f"{len(new_features) - found:.0f} features identified in order {ii+OS:.0f}")
        found = len(new_features)
    return new_features

####################################################################
def thar_auto(dir_name, file_name, OS, X_Order, Y_Order, view):
    #params
    global old_thar
    old_thar = 'thar_last.fits'         #name of last thar with good dispersion function
    global old_thar_features
    old_thar_features = 'thar_last.dat'  #name of features list for last thar
    global line_list
    line_list = 'thar.dat'               #name of full features list for thar lamp
    global FWHM
    FWHM = 4                                        #approximate FWHM of profile
    global threshold
    threshold = 3                               #threshold for automatic search of features, 0.1 for FOX
    global tolerance
    tolerance = 0.05                               #tolerance in A for auto identification of features
    global max_shift
    max_shift = 20                                  #max shift in pix for automatic identification
    thar = []           #full list with thar faetures
    zero_features = []  #list with features for first identification

    #open file with last thar
    with fits.open(old_thar) as hdulist:
        zero = hdulist[0].data.copy()
        zero = np.nan_to_num(zero)

    #open file with new thar
    print(file_name)
    logging.info(file_name)
    with fits.open(file_name) as hdulist:
        prihdr = hdulist[0].header
        spectrum = hdulist[0].data.copy()
        if 'EXPTIME' in prihdr:
            exptime = prihdr['EXPTIME']
        elif 'EXPOSURE' in prihdr:
            exptime = prihdr['EXPOSURE']
        else:
            exptime = 1
        # spectrum = spectrum / exptime
        spectrum = np.nan_to_num(spectrum)
    order = 0

    #read short file with features
    if  os.path.exists(old_thar_features):
        print('File ', old_thar_features, ' found')
        logging.info(f"File {old_thar_features} found")
        with open(old_thar_features, 'r') as f:
            for line in f:
                one_str = line.rsplit('\t')
                zero_features.append([float(one_str[0]),float(one_str[1]), float(one_str[2]), float(one_str[3])])
        f.close()
        print(len(zero_features), "features read for auto identification")
        logging.info(f"{len(zero_features)} features read for auto identification")
    else:
        print('File ', old_thar_features, ' not found')
        logging.error(f"File {old_thar_features} not found")
        return None

    # Read full list of ThAr lines
    if  os.path.exists(line_list):
        print('File ', line_list, ' found')
        logging.info(f"File {line_list} found")
        with open(line_list, 'r') as f:
            for line in f:
                try:
                    thar.append(float(line))
                except:
                    pass
        f.close()
        print(f"{len(thar)} features in line list")
        logging.info(f"{len(thar)} features in line list")

    else:
        print('File ', line_list, ' not found')
        logging.info(f"File {line_list} not found")

    # identifying and centering features from a short list
    new_features = first_ident(spectrum, zero, zero_features, OS)
    # write a new short list with identifications
    write_new_short(spectrum, new_features, prihdr, dir_name)
    print('New short list saved')
    logging.info('New short list saved')
    # search for the first solution based on the short list
    disp_params, new_features, points, rms = dispers_p(new_features, X_Order, Y_Order, view)
    # add new lines
    new_features = add_lines(spectrum, OS, new_features, thar, disp_params, Y_Order)
    # now search for the second solution based on the long list
    disp_params, new_features, points, rms = dispers_p(new_features, X_Order, Y_Order, view)
    # save the results
    write_disp(file_name, disp_params, OS, X_Order, Y_Order, prihdr)
    hdulist[0].data = spectrum
    write_features(file_name, new_features, prihdr)
    hdulist.writeto(file_name, overwrite=True)

    fig = plt.figure(figsize=(16, 5), tight_layout=True)
    ax = fig.add_subplot(111)
    local=np.asarray(new_features)
    O = np.copy(local[:, 0])
    X = np.copy(local[:, 1])
    W = np.copy(local[:, 2])
    residual = solution(disp_params, Y_Order)(X,O)
    residual = W - residual
    rms = round(np.std(residual),5)
    label_data = f"Chebyshev, Xorder={X_Order}, Yorder={Y_Order}, {O.shape[0]} points, RMS={rms}Å  ({np.std((residual/W)*c_const):.2f} km/s)"
    plt.cla()
    ax.plot(O, residual, 'go')
    ax.set_xlim(O[0]-1, O[-1]+1)

    plt.title(label_data)
    plt.ylabel('Angstroms', fontsize=15)
    plt.xlabel('Order', fontsize=15)
    plt.savefig(file_name.replace('.fits', '_disp.pdf'), dpi=350)
    if view:
        plt.show()

    return None
