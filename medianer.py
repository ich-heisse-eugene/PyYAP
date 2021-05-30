import astropy.io.fits as pyfits
import astropy.time as Time
import os
import numpy
import time

##################################################################
def medianer(dir_name, list_name, out_name):
    _data = []
    file_list = []
    _time = []
    JD1 = []
    JD2 = []
    JD3 = []

    ##collect data and time
    with open(list_name, 'r') as f:
        for line in f:
            name = line.strip()
            _data.append(pyfits.getdata(name))
            prihdr = pyfits.getheader(name)
            file_list.append(name.split(os.sep)[-1])
            try:
                _time.append(time.mktime(time.strptime(prihdr['DATE-OBS'], "%Y-%m-%dT%H:%M:%S.%f")))
            except:
                _time.append(time.mktime(time.strptime(prihdr['DATE-OBS']+'.0', "%Y-%m-%dT%H:%M:%S.%f")))
            tm = Time.Time(prihdr['DATE-OBS'])
            JD1.append(tm.jd)
            JD2.append(tm.jd + prihdr['EXPTIME']/2./86400.)
            JD3.append(tm.mjd + prihdr['EXPTIME']/2./86400.)
    f.close()

    ##mean time
    _time = numpy.mean(_time)
    mean_time = time.strftime('%Y-%m-%dT%H:%M:%S.000', time.localtime(_time))
    JD1 = numpy.mean(JD1)
    JD2 = numpy.mean(JD2)
    JD3 = numpy.mean(JD3)
    ##median data
    _data=numpy.asarray(_data)
    out = numpy.median(_data, 0)
    out = numpy.float32(out)

    if out_name == 's_bias.fits':
        RN = []
        for ii in range(0,_data.shape[0]):
            xxx = (numpy.std(_data[int(ii), int(_data.shape[1]/2-10):int(_data.shape[1]/2+10), int(_data.shape[2]/2-10):int(_data.shape[2]/2+10)]))
            RN.append(xxx)
        RN =numpy.asarray(RN)
        ReadNoise = numpy.median(RN)

    mean = numpy.mean(out)
    median = numpy.median(out)
    stdv = numpy.std(out)

    prihdr['DATE-OBS'] = (mean_time, 'mean time of combined files')
    prihdr['JD-OBS']  	= (JD1, 'JD of start of exposure')
    prihdr['JDMIDDLE']	= (JD2, 'JD of middle of exposure')
    prihdr['MJD-OBS'] 	= (JD3, 'MJD of start of exposure')
    prihdr['DATAMAX'] = (numpy.max(out), 'Max pixel value')
    prihdr['DATAMIN'] = (numpy.min(out), 'Min pixel value')
    prihdr['DATAMEAN'] = (mean, 'Mean value')

    if out_name == 's_bias.fits':
        ReadNoise = round(ReadNoise*prihdr['GAIN'],2)
        prihdr['RDNOISE'] = (ReadNoise, 'measured read noise (e)')

    prihdr['HISTORY'] = 'median combined: '+ str(_data.shape[0])

    hdu = pyfits.PrimaryHDU(out, prihdr)
    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto(dir_name.joinpath(out_name), clobber = True)

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
