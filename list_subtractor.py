from astropy.io import fits
import os
import numpy

import warnings
warnings.simplefilter("ignore")
##################################################################
def list_subtractor(list_name, subtrahend_name, stype):
    hdulist = fits.open(subtrahend_name)
    delta = hdulist[0].data.copy()
    prihdr = hdulist[0].header
    hdulist.close()

    with open(list_name, 'r') as f:
        for line in f:
            name = line.strip()
            hdulist = fits.open(name, mode = 'update')
            data = hdulist[0].data.copy()
            prihdr = hdulist[0].header
            hdulist.close()

            if data.shape[0]==delta.shape[0] and data.shape[1]==delta.shape[1]:
                data=numpy.float32(data)-delta
                prihdr['HISTORY'] = stype+' subtracted'
                hdu = fits.PrimaryHDU(data, prihdr)
                hdulist = fits.HDUList([hdu])
                hdulist.writeto(name, overwrite=True)
            else:
                print ("Frame", name, "has wrong size")
    f.close()

    return ("Cleaned")
