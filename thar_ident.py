import math
from astropy.io import fits
import os
import numpy
from sys import argv

from pylab import *
import scipy.ndimage
from scipy.optimize import curve_fit
from scipy import optimize

import matplotlib
import matplotlib.pyplot

import os.path

##import tkinter as Tk
matplotlib.pyplot.ion()
##################################################################
####parameters
base_name = 'thar.dat'
FWHM = 6
OS = 36         # Absolute number of the first (red) order
threshold = 0.001
tolerance = 0.05
X_Order = 5
Y_Order = 5

##################################################################
#define hot keys
key_mark = 'c'      #mark feature and set wavelength
key_del = 'x'       #delete feature
key_up = 'up'       #next order
key_down = 'down'   #prev order
key_fit = 'r'       #fit data
key_fit_o = 'o'     #fit data
key_quit = 'e'      #quit from fit mode
key_auto = 'a'      #add lines auto
key_param = ':'     #enter parameters ([x_order], [y_order], [tolerance] for auto search (A), [OS] order shift, [threshold] for auto search, [FWHM] of features)
key_print = 'ctrl+p'#print params
key_write = 'w'     #write data to files
key_help = '?'
key_grid = 'g'      #grid on


##################################################################
### enter parameter
def param_enter():
    global FWHM
    global OS
    global threshold
    global tolerance
    global X_Order
    global Y_Order
    global fit_params
    global loc_area
    global fit_area
    global features

    param = input("Enter [parameter=X.XX]")
    one_str = param.rsplit('=')
    if one_str[0]== 'FWHM' or one_str[0]== 'fwhm':
        FWHM=float(one_str[1])
        print ('FWHM = ', FWHM)
    elif one_str[0]== 'OS' or one_str[0]== 'shift':
        OS_n=float(one_str[1])
        for i, x in enumerate(features):
            features [i] = (x[0]+(OS_n-OS), x[1], x[2], x[3])
        OS = OS_n
        print ('Order shift = ', OS)
    elif one_str[0]== 'threshold':
        threshold=float(one_str[1])
        print ('Search threshold = ', threshold)
    elif one_str[0]== 'tolerance':
        tolerance=float(one_str[1])
        print ('Wave tolerance = ', tolerance)
    elif one_str[0]== 'X_Order' or one_str[0]== 'x_order':
        X_Order=int(one_str[1])
        print ('X_Order = ', X_Order)
    elif one_str[0]== 'Y_Order' or one_str[0]== 'y_order':
        Y_Order=int(one_str[1])
        print ('Y_Order = ', Y_Order)
    else:
        print('Wrong command. Press shift+? for help')
    fit_params = numpy.zeros([1+X_Order*Y_Order])
    loc_area = int(FWHM*2)
    fit_area = int(FWHM*0.7)

##################################################################
###
def print_help():
    print(
        'c - mark feature and set wavelength', '\n',
        'x - delete feature', '\n',
        'up - next order', '\n',
        'down - prev feature', '\n',
        'r - fit data', '\n',
        'q - quit from fit mode', '\n',
        'a - add lines auto', '\n',
        ': - enter parameters', '\n',
        '[x_order]', '\n',
        '[y_order]', '\n',
        '[tolerance] for auto search (A)', '\n',
        '[OS] order shift', '\n',
        '[threshold] for auto search', '\n',
        '[FWHM] of features', '\n',
        'ctrl+p - print parameters', '\n',
        'w - write data to files', '\n',
        'g - grid on'
            )

##################################################################
###
def print_status():
    print ('FWHM = ', FWHM)
    print ('Order shift = ', OS)
    print ('Search threshold = ', threshold)
    print ('Wave tolerance = ', tolerance)
    print ('X_Order = ', X_Order)
    print ('Y_Order = ', Y_Order)
    print ('Mode = ', mode)
    print ('Status = ', status)

##################################################################
def grid():
    grid = matplotlib.pyplot.figure(2)
    grid_ax = grid.add_subplot(111)
    ii=0
    shift = fit_params[0]
    qqq = numpy.delete(fit_params,0)
    coeff = numpy.reshape(qqq,(-1,Y_Order))
    for ii in range (int(OS), int(OS+spectrum.shape[0])):
        X = numpy.linspace(0, 2047)
        Y = (numpy.polynomial.chebyshev.chebval2d(X,ii,coeff)-shift)/ii
        grid_ax.plot(X,Y)

    local=numpy.asarray(features)
    X=copy(local[:, 1])
    W=copy(local[:, 2])
    grid_ax.plot(X,W, 'go')
    grid_ax.grid()
    matplotlib.pyplot.draw()
##################################################################
### 2D dispersion function fitting
def solution(C): #2d chebyshev polynom calculation
    shift = C[0]
    C = numpy.delete(C,0)
    coeff = numpy.reshape(C,(-1,Y_Order))#7-Yorder
    return lambda X,O:(numpy.polynomial.chebyshev.chebval2d(X,O,coeff)-shift)/O

def dispers_p(local, color):
    local=numpy.asarray(local)
    color=numpy.asarray(color)
    O=copy(local[:, 0])
    X=copy(local[:, 1])
    W=copy(local[:, 2])
    g=numpy.where(color!=0)
    b=numpy.where(color!=1)

    params = numpy.zeros_like(fit_params)
    errorfunction = lambda C: ravel(solution(C)(X,O) - W)
    C, success = optimize.leastsq(errorfunction, params, maxfev=10000)

    resudial = solution(C)(X,O)
    resudial = W-resudial
    rms = round(numpy.std(resudial),5)
    label_data = 'Chebyshev, Xorder='+str(X_Order)+', Yorder='+str(Y_Order)+',' + ' points='+ str(O.shape[0])  + ', RMS=' + str(rms)
    matplotlib.pyplot.cla()
    ax.set_xlim([numpy.min(X), numpy.max(X)])
    ax.plot(X[g], resudial[g], 'go')
    ax.plot(X[b], resudial[b], 'bo')
    matplotlib.pyplot.title(label_data)
    matplotlib.pyplot.draw()

    return (C)

def dispers_o(local,color):
    local=numpy.asarray(local)
    color=numpy.asarray(color)
    O=copy(local[:, 0])
    X=copy(local[:, 1])
    W=copy(local[:, 2])
    g=numpy.where(color!=0)
    b=numpy.where(color!=1)

    params = numpy.zeros_like(fit_params)
    errorfunction = lambda C: ravel(solution(C)(X,O) - W)
    C, success = optimize.leastsq(errorfunction, params, maxfev=10000)

    resudial = solution(C)(X,O)
    resudial = W-resudial
    rms = round(numpy.std(resudial),5)
    label_data = 'Chebyshev, Xorder='+str(X_Order)+', Yorder='+str(Y_Order)+',' + ' points='+ str(O.shape[0])  + ', RMS=' + str(rms)
    matplotlib.pyplot.cla()
    ax.set_xlim([numpy.min(O)-1, numpy.max(O)+1])
    ax.plot(O[g], resudial[g], 'go')
    ax.plot(O[b], resudial[b], 'bo')
    matplotlib.pyplot.title(label_data)
    matplotlib.pyplot.ylabel('Angstroms', fontsize=15)
    matplotlib.pyplot.xlabel('Order', fontsize=15)
    matplotlib.pyplot.draw()

    return (C)

####################################################################
### checker
def WL_checker(WL, o):
    for ii in range (0, len(features)):
            if WL ==features[ii][2] and o+OS==features[ii][0]:
                return (1)
    return (0)

####################################################################
###1D Moffat fitting/center search
def center(x_coo, y_coo, bckgrnd):
    ROI =  copy(spectrum[order, int(x_coo-fit_area):int(x_coo+fit_area)])
    X=numpy.arange(0, ROI.shape[0])
    x0 = int(ROI.shape[0]/2)
    #moffat fitting
    #A - amplitude, B,C - coeff of function (width ...), D - background
    moffat = lambda x, A, B, C, D, x0: A*(1 + ((x-x0)/B)**2)**(-C)+D
    p0 = np.array([y_coo, 3, 3, bckgrnd, x0])
    try:
        popt, pcov = curve_fit(moffat, X, ROI, p0, maxfev=10000)
    except RuntimeError:
        return(0)
    return (popt[4]-fit_area)

####################################################################
###search max element in 1d array
def get_max(x_coo):
    ROI =  copy(spectrum[order, int(x_coo-FWHM):int(x_coo+FWHM)])
    Max = numpy.amax(ROI)
    Min = numpy.amin(ROI)
    for ii in range (0, ROI.shape[0]):
        if ROI[ii]==Max:
            return (x_coo+ii-FWHM, Max, Min)

####################################################################
###redraw spectrum
def draw_spec(spectrum, line):
    graph_data = spectrum[line,:]
    label_data = 'Beam =' + str(line)+ ', '+ 'Order = ' + str(line+OS)+', Status = ' + status

    matplotlib.pyplot.cla()
    ax.plot(graph_data, label=label_data, color= 'black', gid='object')
##    if status=='fit':
##        ax_w = ax.twiny()
##        new_tick_locations = np.array([0.0, 500.0, 1000.0, 1500.0, 2000.0])
##        ax_w.set_xticks(new_tick_locations)
##        ax_w.set_xticklabels(tick_function(new_tick_locations,line+69))
    ax.set_xlim([0, len(graph_data)])
    matplotlib.pyplot.title(label_data)
    if len(features)!=0:
        for ii in range(0,len(features)):
            if features[ii][0]==line+OS:
                x_coo = features[ii][1]
                WL = features[ii][2]
                y_coo = features[ii][3]
                col='g^'
                if color[ii]==0:
                    col='b^'
                ax.plot(x_coo, 0, col)
                ax.annotate(WL, xy=(x_coo, 0), rotation=90,\
                            horizontalalignment='center',\
                            verticalalignment='top',  color='b')
    matplotlib.pyplot.draw()

####################################################################
### mark feature
def mark_feature(x_coo, line):
    x_coo = x_coo.round()
    offset = get_max(x_coo)#search max pixel near click
    print('Max=', offset)
    x_coo = offset[0]
    y_coo = offset[1]
    bckgrnd = offset[2]
    offset = center(x_coo, y_coo, bckgrnd)
    if offset!=0:
        x_coo=round(x_coo+offset,3)
        print('x=',x_coo)
    if status=='unfit':
        WL = input("Enter wavelegnth: ")
    else:
        PWL = solution(fit_params)(x_coo,line+OS)
        loc_thar = (thar - PWL)**2
        nearest = numpy.argmin(loc_thar)
        txt = "Enter wavelegnth["+str(thar[nearest])+"]:"
        WL = input(txt) or thar[nearest]
    try:
        WL = float(WL)
        if WL_checker(WL, order)==0:
            features.append([line+OS, x_coo, WL, y_coo])
            color.append(0)
            ax.plot(x_coo, y_coo, 'bo')
            ax.annotate(WL, xy=(x_coo,y_coo+1),\
                        rotation=90,horizontalalignment='center',  \
                        verticalalignment='baseline',  \
                        color='b')
            matplotlib.pyplot.draw()
    except:
        print('error')
        pass

####################################################################
### delete feature
def del_feature(x_coo, line):
    local=[]
    local_index=[]
    for ii in range(0, len(features)):
        if features[ii][0] == (line+OS):
            local.append(features[ii][1])
            local_index.append(ii)
    local=numpy.array(local)
    local_index=numpy.array(local_index)
    local=(local-x_coo)*(local-x_coo)
    nearest = numpy.argmin(local)
    del features[local_index[nearest]]
    del color[local_index[nearest]]
    draw_spec(spectrum, order)
    matplotlib.pyplot.draw()

####################################################################
### search local maximum in current order
def auto_search(threshold, tolerance):
    row = copy(spectrum[order,0:spectrum.shape[1]])                                     #copy of one row
    row_S = scipy.ndimage.gaussian_filter(row, int(FWHM/3))                     #smooth
    row_f = []                                                                          #array for auto features
    for jj in range (loc_area,row_S.shape[0]-loc_area):                                                                                        #background
        if (row_S[jj]-row_S[jj-1]>0 and row_S[jj+1]-row_S[jj]<0):                       #search for local extremum of smoothed spectra
            offset = get_max(jj)                                                        #search max pixel near click
            x_coo = offset[0]
            y_coo = offset[1]
            bckgrnd = offset[2]
            if ((y_coo-bckgrnd)>threshold):                                             #check threshold
                PWL = solution(fit_params)(x_coo,order+OS)
                loc_thar = (thar - PWL)**2
                nearest = numpy.argmin(loc_thar)
                if math.sqrt(loc_thar[nearest])<tolerance*2:
                    offset = center(x_coo, y_coo, bckgrnd)
                    if offset!=0:
                        x_coo=round(x_coo+offset,3)
                        PWL = solution(fit_params)(x_coo,order+OS)
                        loc_thar = (thar - PWL)**2
                        nearest = numpy.argmin(loc_thar)
                        if math.sqrt(loc_thar[nearest])<tolerance:
                            WL = thar[nearest]
                            if WL_checker(WL, order)==0:
                                features.append([order+OS, x_coo, WL, y_coo])
                                color.append(1)
                                ax.plot(x_coo, y_coo, 'ro')
                                ax.annotate(WL, xy=(x_coo,y_coo+1),\
                                            rotation=90,\
                                            horizontalalignment='center',\
                                            verticalalignment='baseline',  \
                                            color='b')

    matplotlib.pyplot.draw()
####################################################################
### delete feature
def del_point(x_coo, y_coo):
    global fit_params
    local=numpy.asarray(features)
    O=copy(local[:, 0])
    X=copy(local[:, 1])
    W=copy(local[:, 2])

    resudial = solution(fit_params)(X,O)
    resudial = W-resudial

    if mode=='fit':
        scale = (numpy.max(X) - numpy.min(X))/(numpy.max(resudial) - numpy.min(resudial))
        resudial = resudial*scale
        y_coo=y_coo*scale
        distance = sqrt((resudial-y_coo)**2+(X-x_coo)**2)
        x_coo = features[numpy.argmin(distance)][1]
        y_coo = features[numpy.argmin(distance)][2] - solution(fit_params)(x_coo,features[numpy.argmin(distance)][0])
        ax.plot(x_coo, y_coo, 'ro')
        del features[numpy.argmin(distance)]
        del color[numpy.argmin(distance)]

    if mode=='fit_o':
        X=O
        scale = (numpy.max(X) - numpy.min(X))/(numpy.max(resudial) - numpy.min(resudial))
        resudial = resudial*scale
        y_coo=y_coo*scale
        distance = sqrt((resudial-y_coo)**2+(X-x_coo)**2)
        x_coo = features[numpy.argmin(distance)][0]
        y_coo = features[numpy.argmin(distance)][2] - solution(fit_params)(features[numpy.argmin(distance)][1],features[numpy.argmin(distance)][0])
        ax.plot(x_coo, y_coo, 'ro')
        del features[numpy.argmin(distance)]
        del color[numpy.argmin(distance)]

    matplotlib.pyplot.draw()

####################################################################
def write_disp(name, data):
    name = os.path.splitext(name)[0] + "_disp.dat"
    text_file = open(name, "w")
    text_file.write(str(OS)+'\n')
    text_file.write(str(X_Order)+'\n')
    text_file.write(str(Y_Order)+'\n')
    for ii in range(0, len(data)):
        text_file.write(str(data[ii])+'\n')
    text_file.close()
    prihdr['WLSLTN'] = (str(name), 'wavelength solution')

####################################################################
def write_features(features_name, color):
    text_file = open(features_name, "w")
    for ii in range(0, len(features)):
        text_file.write(str(int(features[ii][0]))+'\t'+str(round(features[ii][1],2))+'\t'+str(round(features[ii][2],4))+'\t'+str(round(features[ii][3],1))+'\t'+str(color[ii])+'\n')
    text_file.close()
    prihdr['FEATURES'] = (str(features_name), 'ThAr features for WLSLTN')
##    print(prihdr.ascardlist())

####################################################################
### press event handler
def press(event):
    global order
    global mode
    global fit_params
    global status
    global features_name
    global color
    try:
        if 'alt' in event.key:
            print('Warning! Numlock on!')
    except:
        pass

    if mode=='mark':
        if event is not None and event.key==key_up:
            order = order+1
            if order>spectrum.shape[0]-1:
                order=0
            draw_spec(spectrum, order)

        if event is not None and event.key==key_down:
            order = order-1
            if order<0:
                order=spectrum.shape[0]-1
            draw_spec(spectrum, order)

        if event is not None and event.key==key_mark:
            mark_feature(event.xdata, order)

        if event is not None and event.key==key_auto:
            if status=='fit':
                auto_search(threshold, tolerance)

        if event is not None and event.key==key_del:
            del_feature(event.xdata, order)

        if event is not None and event.key==key_param:
            param_enter()

        if event is not None and event.key==key_print:
            print_status()

        if event is not None and event.key==key_help:
            print_help()

        if event is not None and event.key==key_write:
            try:
                fit_params = dispers_p(features, color)
                write_disp(features_name, fit_params)
                mode='fit'
            except:
                print('unable fit')
                pass
            write_features(features_name,color)
            hdulist[0].data = spectrum
##            print(prihdr.ascardlist())
            hdulist.writeto(fits_name, overwrite=True)

        if event is not None and event.key==key_fit:
##            try:
            fit_params = dispers_p(features, color)
            mode='fit'
##            except:
##                print('unable fit')
##                pass

        if event is not None and event.key==key_fit_o:
            try:
                fit_params = dispers_o(features, color)
                mode='fit_o'
            except:
                print('unable fit')
                pass

    elif mode=='fit' or mode=='fit_o':
        if event is not None and event.key==key_quit:
            mode='mark'
            status = 'fit'
            draw_spec(spectrum, order)

        if event is not None and event.key==key_fit:
            fit_params = dispers_p(features, color)

        if event is not None and event.key==key_fit_o:
            fit_params = dispers_o(features,color)

        if event is not None and event.key==key_del:
            del_point(event.xdata, event.ydata)

        if event is not None and event.key==key_param:
            param_enter()

        if event is not None and event.key==key_print:
            print_status()

        if event is not None and event.key==key_write:
            fit_params = dispers_p(features,color)
##            write_disp(features_name, fit_params)
            write_features(features_name, color)
            hdulist[0].data = spectrum
##            print(prihdr.ascardlist())
            hdulist.writeto(fits_name, overwrite=True)

        if event is not None and event.key==key_grid:
            grid()

##        print ("Quit")
##        matplotlib.pyplot.close('all')
##        os._exit(0)

####################################################################
def thar_manual(file_name):

    #read file with spectrum
    global fits_name
    global hdulist
    global spectrum
    global order
    global features_name
    global prihdr
    global color

    fits_name = file_name
    features_name = os.path.splitext(file_name)[0] + '.dat'

    hdulist = fits.open(file_name)
    prihdr = hdulist[0].header
    spectrum = hdulist[0].data.copy()
    expthar = prihdr['EXPTIME']
    print(f"Texp = {expthar} s")
    spectrum = spectrum / expthar
    hdulist.close()
    order = 0

    #read features file
    if  os.path.exists(features_name):
        print('File ', features_name, ' found')
        with open(features_name, 'r') as f:
            for line in f:
                one_str = line.rsplit('\t')
                if len(one_str)==5:
                    features.append([float(one_str[0]),float(one_str[1]), float(one_str[2]), float(one_str[3])])
                    color.append(float(one_str[4]))
                else:
                    features.append([float(one_str[0]),float(one_str[1]), float(one_str[2]), float(one_str[3])])
                    color.append(1)
        f.close()
        print(len(features), "features read from file")

    else:
        print('File ', features_name, ' not found')

    #read thar file
    if  os.path.exists(base_name):
        print('File ', base_name, ' found')
        with open(base_name, 'r') as f:
            for line in f:
                try:
                    thar.append(float(line))
                except:
                    pass
        f.close()
        print(len(thar), "features in line list")

    else:
        print('File ', base_name, ' not found')

    cid = fig.canvas.mpl_connect('key_press_event', press)
    draw_spec(spectrum, order)

    return (None)

file_thar = argv[1]

fig = matplotlib.pyplot.figure(1, figsize=(16, 5))
ax = fig.add_subplot(111)

features=[]
color=[]
thar = []
fit_params = numpy.zeros([1+X_Order*Y_Order])
loc_area = int(FWHM*2)
fit_area = int(FWHM*0.7)
mode = 'mark'
status = 'unfit'
features_name=''

thar_manual(file_thar)
matplotlib.pyplot.ioff()
matplotlib.pyplot.show()
