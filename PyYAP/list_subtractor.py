import astropy.io.fits as pyfits
import os
import numpy

import warnings
warnings.simplefilter("ignore")
##################################################################
def list_subtractor(dir_name, list_name, subtrahend_name):
    hdulist = pyfits.open(subtrahend_name)
    delta = hdulist[0].data.copy()
    prihdr = hdulist[0].header
    hdulist.close()

    with open(list_name, 'r') as f:
        for line in f:
            name = line.strip()
            hdulist = pyfits.open(name, mode = 'update')
            data = hdulist[0].data.copy()
            prihdr = hdulist[0].header
            hdulist.close()

            if data.shape[0]==delta.shape[0] and data.shape[1]==delta.shape[1]:
                data=numpy.float32(data)-delta
                prihdr['HISTORY'] = 'bias subtracted'
                hdu = pyfits.PrimaryHDU(data, prihdr)
                hdulist = pyfits.HDUList([hdu])
                hdulist.writeto(name, clobber=True)
            else:
                print ("Frame", name, "has wrong size")
    f.close()

    return ("Cleaned")

##dir_name ='D:/Observations/Melnikov/zTest/'
##list_name ='thar_list.txt'
##list_subtractor(dir_name, list_name, 's_bias.fits')
