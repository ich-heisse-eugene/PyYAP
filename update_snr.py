from astropy.io import fits
import numpy as np
import logging
from get_sp_resolv import read_multispec

def update_snr(input_file, queue):
    file_err = input_file.replace('_ec_', '_err_')
    wl, err = read_multispec(file_err)
    norders = np.shape(wl)[0]
    npix = np.shape(wl)[1]
    pix = np.linspace(0, npix, 5, dtype=int)
    with fits.open(input_file, mode = 'update') as hdulist:
        prihdr = hdulist[0].header
        # Update information about sp. resolution
        refspec = prihdr['REFSPEC'].replace('_disp.txt', '_WCS.fits')
        with fits.open(refspec) as hduref:
            hdrref = hduref[0].header
            R = hdrref['R']
            prihdr.set('RESOL', int(R), 'Median spectral resolution of data')

        #end
        for i in range(norders):
            key_n = 'SNR'+str(i)
            key_v = ' '
            key_c = 'SNR at '
            for x in pix[1:-1]:
                w0 = np.mean(wl[i, x-10:x+10])
                snr0 = np.mean(err[i, x-10:x+10])
                if np.isnan(snr0):
                    print(f"NaN values of SNR in a file {input_file}")
                    queue.put((logging.INFO, f"NaN values of SNR in a file {input_file}"))
                    key_v = key_v + 'NaN'+' '
                else:
                    key_v = key_v + str(int(snr0))+' '
                key_c = key_c + str(int(w0))+'A  '
            prihdr.set(key_n, key_v, key_c)
        hdulist[0].header = prihdr
    return None

