from datetime import datetime, timedelta
import astropy.io.fits as pyfits
import numpy as np

import logging

def thar_manager(dir_name, obj_name, thar_list):

##    try:
    hdulist = pyfits.open(obj_name)
    data = hdulist[0].data.copy()
    prihdr = hdulist[0].header
    hdulist.close()
    obs_time = datetime.strptime(prihdr['DATE-OBS'], "%Y-%m-%dT%H:%M:%S.%f") # Normally 'DATE-OBS'
    exptime = timedelta(seconds = int(prihdr['EXPOSURE']))
    middle_time = obs_time + exptime/2.
    print(f"File: {obj_name}, Middle: {datetime.strftime(middle_time, '%Y-%m-%d %H:%M:%S')[:-3]}")
    logging.info(f"File: {obj_name}, Middle: {datetime.strftime(middle_time, '%Y-%m-%d %H:%M:%S')[:-3]}")

    thar_name = []
    thar_time = []

    with open(thar_list, 'r') as f:
        for line in f:
            name = line.strip()
            try:
                prihdr = pyfits.getheader(name)
                thar_name.append(name)
                try:
                    diff_time = datetime.strptime(prihdr['DATE-OBS'], "%Y-%m-%dT%H:%M:%S.%f") - middle_time
                    diff_time = abs(diff_time.total_seconds())
                    thar_time.append(diff_time)
                except:
                    diff_time = datetime.strptime(prihdr['DATE-OBS'], "%Y-%m-%dT%H:%M:%S.%f") - middle_time
                    diff_time = abs(diff_time.total_seconds())
                    thar_time.append(diff_time)
            except IOError:
                print ("Can't open file:", line)
    f.close()
    thar_time = np.asarray(thar_time)
    return (thar_name[np.argmin(thar_time)], thar_time[np.argmin(thar_time)])
