import astropy.io.fits as pyfits
import os
import numpy as np
import scipy.ndimage
from scipy.optimize import curve_fit
from scipy import optimize
import matplotlib
import matplotlib.pyplot
import os.path
import logging

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
        print(str(int(features[ii][0]))+'\t'+str(round(features[ii][1],2))+'\t'+str(round(features[ii][2],4))+'\t'+str(round(features[ii][3],1)), file=text_file)
    text_file.close()
    prihdr['FEATURES'] = (str(f_name.split(os.sep)[-1]), 'ThAr features for WLSLTN')
    return(prihdr)

####################################################################
def write_new_short(spectrum, features, prihdr, dir_name):
    text_file = open(dir_name.joinpath('thar_new.dat'), "w")
    for ii in range(0, len(features)):
        print(str(int(features[ii][0]))+'\t'+str(round(features[ii][1],2))+'\t'+str(round(features[ii][2],4))+'\t'+str(round(features[ii][3],1)), file=text_file)
    text_file.close()
    hdu = pyfits.PrimaryHDU(spectrum)
    hdu.header = prihdr
    hdu.writeto(dir_name.joinpath('s_thar_new.fits'), overwrite=True)

####################################################################
### checker
def WL_checker(WL, o, OS, features):
    for ii in range (0, len(features)):
            if WL ==features[ii][2] and o+OS==features[ii][0]:
                return (1)
    return (0)

####################################################################
###1D Moffat fitting/center search
def center(s_order, x_coo, y_coo, bckgrnd):
    fit_area = int(FWHM*0.7)
    ROI =  np.copy(s_order[int(x_coo-fit_area):int(x_coo+fit_area)])
    X=np.arange(0, ROI.shape[0])
    x0 = int(ROI.shape[0]/2)
    #moffat fitting
    #A - amplitude, B,C - coeff of function (width ...), D - background
    moffat = lambda x, A, B, C, D, x0: A*(1 + ((x-x0)/B)**2)**(-C)+D
    p0 = np.array([y_coo, 3, 3, bckgrnd, x0])
    try:
        popt, pcov = curve_fit(moffat, X, ROI, p0, maxfev=10000)
    except:
        return(0)
    return (popt[4]-fit_area)

####################################################################
###search max element in 1d array
def get_max(s_order, x_coo):
    width = 1.5*FWHM
    ROI =  np.copy(s_order[int(x_coo-width):int(x_coo+width)])
    Max = np.amax(ROI)
    Min = np.amin(ROI)
    ii = np.unravel_index(ROI.argmax(), ROI.shape)[0]
    return (x_coo+ii-width, Max, Min)

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
    print(f"RMS(km/s): {round(np.std((residual/W)*300000.0),5):.2f}")
    logging.info(f"Lines used: {points}")
    logging.info(f"RMS(Angstrom): {rms:.4f}")
    logging.info(f"RMS(km/s): {round(np.std((residual/W)*300000.0),5):.2f}")

    if view:
        fig = matplotlib.pyplot.figure(1, figsize=(16, 5))
        ax = fig.add_subplot(111)
        points = O.shape[0]
##        residual = (residual/W)*300000.0
        rms = round(np.std(residual),5)
        label_data = 'Chebyshev, Xorder='+str(X_Order)+', Yorder='+str(Y_Order)+',' + ' points='+ str(points)  + ', RMS=' + str(rms)
        matplotlib.pyplot.cla()
        ax.set_xlim([np.min(X), np.max(X)])
        ax.plot(X, residual, 'go')
        matplotlib.pyplot.title(label_data)
        matplotlib.pyplot.ylabel('Angstrom', fontsize=15)
        matplotlib.pyplot.xlabel('Pixel', fontsize=15)
        matplotlib.pyplot.show()

    return (C, features, points, rms)

####################################################################
### auto search lines
def add_lines(spectrum, OS, new_features, thar, disp_params, Y_Order):
    for ii in range(0, spectrum.shape[0]):
        old_len = len(new_features)
        row = spectrum[ii,:]                                                                #copy one order
        row_S = scipy.ndimage.filters.gaussian_filter(row, int(FWHM/3))                     #smooth
        row_f = []                                                                          #array for auto features
##        loc_area = int(FWHM*2)
        for jj in range (int(FWHM*2),row_S.shape[0]-int(FWHM*2)):
            if (row_S[jj]-row_S[jj-1]>0 and row_S[jj+1]-row_S[jj]<0):                       #search for local extremum of smoothed spectra
                offset = get_max(row, jj)                                                   #search max pixel near click
                x_coo = offset[0]
                y_coo = offset[1]
                bckgrnd = offset[2]
                if ((y_coo-bckgrnd)>threshold):                                             #check threshold
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
def search_shift(s_order, z_order,a):
    shift=[]
    for ii in range(-max_shift, max_shift):
        s_shift=ii
        s_data = s_order
        s_x = np.arange(0, len(s_data), 1)
        s_x = s_x+s_shift

        z_data = z_order
        z_x = np.arange(0, len(z_data), 1)

        #check and fix arrays length
        if (len(z_data)-len(s_data))>0:
            add = np.zeros(len(z_data)-len(s_data))
            s_data = np.append(s_data, add)

        if (len(z_data)-len(s_data))<0:
            add = np.zeros(len(s_data)-len(z_data))
            z_data = np.append(z_data, add)

        #add shift
        if s_shift>0:
            add = np.zeros(s_shift)
            z_data = np.append(z_data, add)
            s_data = np.append(add, s_data)
        else:
            add = np.zeros(-s_shift)
            s_data = np.append(s_data, add)
            z_data = np.append(add, z_data)

        #trim zeros
        index = np.arange(0, len(add), 1, int)
        z_data = np.delete(z_data, index)
        s_data = np.delete(s_data, index)
        index = np.arange(len(z_data)-len(add), len(z_data), 1, int)
        z_data = np.delete(z_data, index)
        s_data = np.delete(s_data, index)
        delta = np.sum((z_data-s_data)*(z_data-s_data))
        shift.append(delta)

    s_shift = np.argmin(shift)-max_shift
    return(s_shift)

####################################################################
###reidentify features from short list
def reidentify_features(s_order, OS, zero_features, line, shift, new_features):
    for ii in range(0,len(zero_features)):
        if zero_features[ii][0]==line+OS:
            WL = zero_features[ii][2]
            x_coo = zero_features[ii][1]
            x_coo_new = x_coo-shift

            try:
                x_coo_new = x_coo_new.round()
                offset = get_max(s_order, x_coo_new)     #search max pixel near click
                x_coo_new = offset[0]
                y_coo_new= offset[1]
                bckgrnd = offset[2]
                offset = center(s_order, x_coo_new, y_coo_new, bckgrnd)
                if offset!=0:
                    x_coo_new=round(x_coo_new+offset,3)
                    new_features.append([line+OS, x_coo_new, WL, y_coo_new])
            except:
                pass
    return(new_features)

######################################################################
def first_ident(spectrum, zero, zero_features, OS):
    a=0
    new_features=[]
    #for every order
    founded=0
    for ii in range(0, min(zero.shape[0], spectrum.shape[0])):#zero.shape[0]):
        s_order = spectrum[ii,:]                #copy one order
        z_order = zero[ii,:]
##        shift = search_shift(s_order, z_order)  #get shift
        if ii >53:
            a=1
        shift = search_shift(s_order, z_order,a)  #get shift
        new_features = reidentify_features(s_order, OS, zero_features, ii, shift, new_features)
        print(f"shift= {str(shift)}\t'{len(new_features) - founded:.0f} features identified in order {ii+OS:.0f}")
        logging.info(f"shift= {str(shift)}\t'{len(new_features) - founded:.0f} features identified in order {ii+OS:.0f}")
        founded = len(new_features)
    return(new_features)

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
    FWHM = 6                                        #approximate FWHM of profile
    # global OS
    # OS = 51                                         #absolute number of first (red) order
    global threshold
    threshold = 0.05                               #threshold for automatic search of features, 0.1 for FOX
    global tolerance
    tolerance = 0.05                               #tolerance in A for auto identification of features
    # global X_Order
    # X_Order = 6                                     #X (echelle dispersion) order of global 2D polynome
    # global Y_Order
    # Y_Order = 6                                     #Y (cross dispersion) order of global 2D polynome
    global max_shift
    max_shift = 20                                  #max shift in pix for automatic identification
    thar = []           #full list with thar faetures
    zero_features = []  #list with features for first identification

    #open file with last thar
    hdulist = pyfits.open(old_thar)
    zero = hdulist[0].data.copy()
    hdulist.close()
    zero = np.nan_to_num(zero)

    #open file with new thar
    print(file_name)
    logging.info(file_name)
    hdulist = pyfits.open(file_name) ## dir_name+file_name (for Windows???)
    spectrum = hdulist[0].data.copy()
    prihdr = hdulist[0].header
    hdulist.close()
    spectrum =np.nan_to_num(spectrum)
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
        return (None)

    #read full thar lines list
    if  os.path.exists(line_list):
        print('File ', line_list, ' founded')
        logging.info(f"File {line_list} founded")
        with open(line_list, 'r') as f:
            for line in f:
                try:
                    thar.append(float(line))
                except:
                    pass
        f.close()
        print(len(thar), "features in line list")
        logging.info(f"{len(thar)} features in line list")

    else:
        print('File ', line_list, ' not found')
        logging.info(f"File {line_list} not found")

    new_features=first_ident(spectrum, zero, zero_features, OS)                 #identify and centering features from short list
    write_new_short(spectrum, new_features, prihdr, dir_name)                         #write new short list for identification
    print('New short list saved')
    logging.info('New short list saved')
    disp_params, new_features, points, rms = dispers_p (new_features, X_Order, Y_Order, view)       #search first solution for short list
    new_features = add_lines(spectrum, OS, new_features, thar, disp_params, Y_Order)                #add lines
    disp_params, new_features, points, rms = dispers_p (new_features, X_Order, Y_Order, view)       #search solution for long list
    write_disp(file_name, disp_params, OS, X_Order, Y_Order, prihdr)
    hdulist[0].data = spectrum
    write_features(file_name, new_features, prihdr)
    hdulist.writeto(file_name, overwrite=True)

    fig = matplotlib.pyplot.figure(figsize=(16, 5), tight_layout=True)
    ax = fig.add_subplot(111)
    local=np.asarray(new_features)
    O = np.copy(local[:, 0])
    X = np.copy(local[:, 1])
    W = np.copy(local[:, 2])
    residual = solution(disp_params, Y_Order)(X,O)
    residual = W-residual
    rms = round(np.std(residual),5)
    label_data = 'Chebyshev, Xorder='+str(X_Order)+', Yorder='+str(Y_Order)+',' + ' points='+ str(O.shape[0])  + ', RMS=' + str(rms) + r'\AA' + f" ({round(np.std((residual/W)*300000.0),5):.2f} km/s)"
    matplotlib.pyplot.cla()
    ax.plot(O, residual, 'go')
    matplotlib.pyplot.title(label_data)
    matplotlib.pyplot.ylabel('Angstroms', fontsize=15)
    matplotlib.pyplot.xlabel('Order', fontsize=15)
    matplotlib.pyplot.savefig(file_name.replace('.fits', '_disp.pdf'), dpi=350)
    if view:
        matplotlib.pyplot.show()

    return (None)
