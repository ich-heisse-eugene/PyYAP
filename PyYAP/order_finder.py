import astropy.io.fits as pyfits
import os
import numpy

from scipy.ndimage.filters import gaussian_filter as GF
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

import matplotlib
import matplotlib.pyplot

"""
Modified by Eugene Semenko
"""

def fit_peak(slic, peaks, cpeak, npeaks):
    if cpeak != npeaks-1:
        half_width = (peaks[cpeak+1,0] - peaks[cpeak,0])/2
    else:
        half_width = (peaks[cpeak,0] - peaks[cpeak-1,0])/2

    half_width = half_width+1 #small correction

    low = int(half_width)
    high = int(half_width)

    if peaks[cpeak,0]-low < 0:
        low = peaks[cpeak,0] - 1
    if peaks[cpeak,0]+high > slic.shape[0]:
        high = slic.shape[0]-peaks[cpeak,0]-1

    local_area = slic[(int(peaks[cpeak,0]-low)) : (int(peaks[cpeak,0]+high))]
    x0=low
    X=numpy.arange(0, local_area.shape[0])

    #moffat fitting
    #A - amplitude, B,C - coeff of function (width ...), D - background
    moffat = lambda x, A, B, C, D, x0: A*(1 + ((x-x0)/B)**2)**(-C)+D
    p0 = numpy.array([peaks[cpeak,1], 3, 3, 0, x0])
    try:
        popt, pcov = curve_fit(moffat, X, local_area, p0, maxfev=10000)
    except RuntimeError:
        pass
    FWHM = 2*popt[1]*numpy.sqrt(2**(1/(popt[2]))-1)
    return local_area, popt, FWHM, low

##################################################################
def order_finder(dir_name, file_name, threshold, slice_width, view):
    orders = []

    #read file
    hdulist = pyfits.open(dir_name+'/'+file_name)
    arr = hdulist[0].data.copy()

    #slice near center column
    center_slice = arr[:,int(arr.shape[1]/2-slice_width):int(arr.shape[1]/2+slice_width)]
    #median
    center_slice = numpy.median(center_slice, 1)
    # #smooth
    # center_slice_S = GF(center_slice, 3)

    #search local maximum
    peaks_coord,_ = find_peaks(center_slice, height=50)
    peaks = numpy.vstack((peaks_coord, center_slice[peaks_coord], numpy.zeros(len(peaks_coord)))).transpose()

    #plot slice and mark founded orders
    if view:
        pix = numpy.arange(center_slice.shape[0])
        fig_center_slice = matplotlib.pyplot.figure(2)
        ax_center_slice = fig_center_slice.add_subplot(111)
        matplotlib.pyplot.cla()
        ax_center_slice.plot(pix, center_slice, 'r')
        ax_center_slice.plot(peaks[:,0], peaks[:,1], 'go')
        ax_center_slice.set_xlim(0, numpy.shape(arr)[0])

    #centering and fit every order
    for ii in range (0, peaks.shape[0]):
        #copy subarray
        local_area, fitpar, FWHM, low  = fit_peak(center_slice, peaks, ii, peaks.shape[0])
        peaks[ii,0] = peaks[ii,0] - low + fitpar[4]
        peaks[ii,1] = fitpar[0]
        peaks[ii,2] = FWHM
        #print (FWHM)

        #plot slice
        if view:
            X=numpy.arange(0, local_area.shape[0], 0.1)
            moffat = lambda x, A, B, C, D, x0: A*(1 + ((x-x0)/B)**2)**(-C)+D
            M_flux=moffat(X, *fitpar)
##            ax_center_slice.plot(peaks[ii,0]-low+X, M_flux)
            ax_center_slice.errorbar(peaks[ii,0], 0, xerr = FWHM, yerr = 0, marker='.')

    if view:
        matplotlib.pyplot.show()

    print (peaks.shape[0], "orders found")
    FWHM = numpy.median(peaks[:, 2])
    print ("Median FWHM: %.2f" %FWHM)
    numpy.savetxt(dir_name +'/temp/orders.txt', peaks)
    print ("Data saved to " + dir_name +'/temp/orders.txt')
    print()
    return dir_name +'/temp/orders.txt'

############################test######################################
# threshold = 60
# slice_width = 5
# view = True
# dir_name = '/home/eugene/pipelines/PyYAP/20200112'
# file_name = 'temp/s_flat.fits'
#
# order_finder(dir_name, file_name, threshold, slice_width, view)
##
