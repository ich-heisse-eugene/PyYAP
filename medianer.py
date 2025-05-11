from astropy.io import fits
import astropy.time as Time
from astropy import units as u
import os
import numpy as np
import time

##################################################################
def medianer(dir_name, list_name, out_name):
    _data = []
    file_list = []
    JD1 = []
    JD2 = []
    JD3 = []

    ##collect data and time
    with open(list_name, 'r') as f:
        for line in f:
            name = line.strip()
            _data.append(fits.getdata(name))
            prihdr = fits.getheader(name)
            file_list.append(name.split(os.sep)[-1])
            tm = Time.Time(prihdr['DATE-OBS'], scale='utc', format='fits')
            print(line, prihdr['DATE-OBS'], tm)
            JD1.append(tm.jd)
            JD2.append(tm.jd + prihdr['EXPTIME']/2./86400.)
            JD3.append(tm.mjd + prihdr['EXPTIME']/2./86400.)

    ##mean time
    JD1 = np.mean(JD1)
    JD2 = np.mean(JD2)
    JD3 = np.mean(JD3)
    mean_time = Time.Time(JD2, format='jd').fits
    ##median data
    _data = np.asarray(_data)
    if str(out_name).find('ordim') != -1 or _data.shape[0] % 2 == 0:
        out = np.mean(_data, 0)
    else:
        out = np.median(_data, 0)
    out = np.float32(out)

    if out_name == 's_bias.fits':
        RN = []
        for ii in range(0,_data.shape[0]):
            area = _data[int(ii), _data.shape[1]//2-10:_data.shape[1]//2+10, _data.shape[2]//2-10:_data.shape[2]//2+10].flatten()
            for cnt in range(5):
                area_std = np.std(area)
                filt = np.where((area >= np.mean(area)-area_std) & (area <= np.mean(area)+area_std))
                area = area[filt]
            xxx = np.std(area)
            RN.append(xxx)
        RN =np.asarray(RN)
        ReadNoise = np.median(RN)

    mean = np.mean(out)
    median = np.median(out)
    stdv = np.std(out)

    print(f"Mean time of the series: {mean_time}")
    prihdr['DATE-OBS'] = (mean_time, 'mean time of combined files')
    prihdr['JD-OBS']  	= (JD1, 'JD of start of exposure')
    prihdr['JDMIDDLE']	= (JD2, 'JD of middle of exposure')
    prihdr['MJD-OBS'] 	= (JD3, 'MJD of start of exposure')
    prihdr['DATAMAX'] = (np.max(out), 'Max pixel value')
    prihdr['DATAMIN'] = (np.min(out), 'Min pixel value')
    prihdr['DATAMEAN'] = (mean, 'Mean value')

    if out_name == 's_bias.fits':
        ReadNoise = round(ReadNoise*prihdr['GAIN'],2)
        prihdr['RDNOISE'] = (ReadNoise, 'measured read noise (e)')

    prihdr['HISTORY'] = 'median combined: '+ str(_data.shape[0])

    hdu = fits.PrimaryHDU(out, prihdr)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(os.path.join(dir_name, out_name), overwrite = True)

    if out_name == 's_bias.fits' or out_name == 's_flat.fits':
        with open(list_name, 'r') as f:
            for line in f:
                name = line.strip()
                try:
                    os.remove(name)
                except IOError:
                    print(f"Can't delete file: {name}")
        f.close()
        os.remove(list_name)

    return(mean, median, stdv)
