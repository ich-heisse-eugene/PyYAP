import astropy.io.fits as pyfits
import os
import numpy as np
from PyAstronomy.pyasl import helcorr
from astropy import time, coordinates as coord, units as u

import logging

import matplotlib
import matplotlib.pyplot
from scipy.interpolate import interp1d
# adds dispersion solution to the fits header
# and linearizes spectrum

c = 299792.458 # km/s

def cheb_sol(C, y_order): #2d chebyshev polynom calculation
    shift = C[0]
    C = np.delete(C,0)
    coeff = np.reshape(C,(-1,y_order))#7-Yorder
    return lambda X,O:(np.polynomial.chebyshev.chebval2d(X,O,coeff)-shift)/O

def disp_add(fits_name, thar_name, view):
    bcr = -999

    if view:
        disp = matplotlib.pyplot.figure(1)
        ax = matplotlib.pyplot.gca()

    #read fits
    hdulist = pyfits.open(fits_name)
    spectrum = hdulist[0].data.copy()
    prihdr = hdulist[0].header
    hdulist.close()

    if 'IMAGETYP' in prihdr:
        if prihdr['IMAGETYP'] == 'OBJ':
            # Compute barycentric correction
            if 'BARYCORR' not in prihdr:
                if 'OBSGEO-B' in prihdr and 'OBSGEO-L' in prihdr and 'OBSGEO-H' in prihdr \
                              and 'RA' in prihdr and 'DEC' in prihdr and 'EPOCH' in prihdr:
                    print("Compute barycentric correction...")
                    logging.info("Compute barycentric correction...")
                    dateobs = prihdr['DATE-OBS']
                    ra = prihdr['RA']
                    dec = prihdr['DEC']
                    epoch = prihdr['EPOCH']
                    obslat = prihdr['OBSGEO-B']
                    obslon = prihdr['OBSGEO-L']
                    obsalt = prihdr['OBSGEO-H']
                    star = coord.SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame='icrs')
                    observat = coord.EarthLocation.from_geodetic(obslon, obslat, obsalt * u.m)
                    dateobs = np.char.replace(dateobs, 'T', ' ')
                    dateobs = time.Time(dateobs, scale='utc', location=observat)
                    ltt_bary = dateobs.light_travel_time(star)
                    bjd = dateobs.jd + ltt_bary.value
                    bcr, hjd = helcorr(observat.lon.degree, observat.lat.degree, obsalt, star.ra.degree, star.dec.degree, dateobs.jd)
                    prihdr.set('HJD', hjd, 'Heliocentric JD')
                    prihdr.set('BJD', bjd, 'Barycentric JD')
                    prihdr.set('BARYCORR', bcr, 'Barycentric correction')
                else:
                    print("No information about the observatory in the header. Wavelength remains uncorrected")
                    logging.warning("No information about the observatory in the header. Wavelength remains uncorrected")
            else:
                bcr = prihdr['BARYCORR']

    #read solution
    sol_name = os.path.splitext(thar_name)[0] + '_disp.txt'
    print('Disp solution:', sol_name)
    solution = np.genfromtxt(sol_name)
    order_shift = int(solution[0])
    x_order = int(solution[1])
    y_order = int(solution[2])
    solution = np.delete(solution, [0, 1, 2])

    WAT2 = ''
    add_string = ''

    X = np.arange(0, spectrum.shape[1], 1)
    for ii in range(0, spectrum.shape[0]):
        O = np.empty_like(X)
        O.fill(ii+order_shift)
        WL = cheb_sol(solution, y_order)(X,O)
        if bcr != -999:
            WL = WL * np.sqrt((1. + bcr/c)/(1. - bcr/c))  # Doppler correction
        w1=np.min(WL)
        w2=np.max(WL)
        nw=len(X)
        dw=(w2-w1)/(nw-1) #-1
##        print ('Blue end:', w1, 'Red end:', w2, 'A/pix:', dw)
        X_lin=np.linspace(w1,w2, nw)
        graph_data = spectrum[ii,:]
        f2 = interp1d(WL, graph_data, kind='slinear', bounds_error = False)
        Y_lin=f2(X_lin)
        spectrum[ii,:] = Y_lin
        APNUM = (prihdr['APNUM'+str(ii+1)])
        del prihdr['APNUM'+str(ii+1)]
        APNUM = APNUM.split()
        add_string = (" spec%i = \"%i %i %i %.10f %0.15f %i %i. %.2f %.2f\"" % (spectrum.shape[0]-ii, spectrum.shape[0]-ii, ii+order_shift, 0, w1, dw, nw, 0, float(APNUM[2]), float(APNUM[3])))
        WAT2 = add_string + WAT2

        if view:
            ax.plot(X_lin, Y_lin, 'r')

    del prihdr['CTYPE1']
    del prihdr['CTYPE2']
    del prihdr['CRPIX1']
    del prihdr['WAT0_001']
    del prihdr['WAT1_001']
    del prihdr['WAT2_001']

    prihdr['CTYPE1'] = 'MULTISPE'
    prihdr['CTYPE2'] = 'MULTISPE'
    prihdr['CRPIX1'] = 0
    prihdr['WAT0_001'] = 'system=multispec'
    prihdr['WAT1_001'] = 'wtype=multispec label=Wavelength units=angstroms'
    WAT2 = 'wtype=multispec' + WAT2
    for ii in range(1, int(2+len(WAT2)/68)):
        keyword = str("WAT2_%03i" % ii)
        prihdr[keyword] =  str(WAT2[0:68])
        WAT2 = WAT2[68:]
    prihdr['REFSPEC'] = (str(sol_name.split(os.sep)[-1]), 'WL reference spectrum')
    if bcr != -999:
        prihdr['HISTORY'] = 'Applied barycentric correction BARYCORR'

    hdulist[0].data = np.flip(spectrum, axis=0)
    new_name = os.path.splitext(fits_name)[0]  + "_WCS.fits"
    hdulist.writeto(new_name, overwrite=True)

    if view:
        matplotlib.pyplot.show()
    return new_name
