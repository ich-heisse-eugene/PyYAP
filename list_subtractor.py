from astropy.io import fits
import os
import numpy
import multiprocessing as mp

import warnings
warnings.simplefilter("ignore")
##################################################################
def process(name, delta, stype):
    with fits.open(name, mode = 'update') as hdulist:
        data = hdulist[0].data.copy()
        prihdr = hdulist[0].header

    if data.shape[0] == delta.shape[0] and data.shape[1] == delta.shape[1]:
        data = numpy.float32(data) - delta
        prihdr['HISTORY'] = stype+' subtracted'
        hdu = fits.PrimaryHDU(data, prihdr)
        hdulist = fits.HDUList([hdu])
        hdulist.writeto(name, overwrite=True)
    else:
        print(f"Image {name} has wrong size")
    return None


def list_subtractor(list_name, subtracted_name, stype, conf):
    with fits.open(subtracted_name) as hdulist:
        delta = hdulist[0].data.copy()
        prihdr = hdulist[0].header

    with open(list_name, 'r') as f:
        nCPUs = os.cpu_count()
        if 'threading' in conf and eval(conf['threading']) and nCPUs > 2:
            proc_args = [(line.strip(), delta, stype) for line in f]
            with mp.Pool(processes=nCPUs) as pool:
                pool.starmap_async(process, proc_args, chunksize=nCPUs)
                pool.close()
                pool.join()
        else:
            for line in f:
                process(line.strip(), delta, stype)
    return("cleaned")
