from datetime import datetime, timedelta
from astropy.io import fits
import numpy as np

import logging

def thar_manager(obj_name, thar_list, queue):

    with fits.open(obj_name) as hdulist:
        data = hdulist[0].data.copy()
        prihdr = hdulist[0].header
    obs_time = datetime.strptime(prihdr['DATE-OBS'], "%Y-%m-%dT%H:%M:%S.%f")
    middle_time = obs_time
    print(f"File: {obj_name}, Middle: {datetime.strftime(middle_time, '%Y-%m-%d %H:%M:%S')[:-3]}")
    queue.put((logging.INFO, f"File: {obj_name}, Middle: {datetime.strftime(middle_time, '%Y-%m-%d %H:%M:%S')[:-3]}"))

    thar_name = []
    thar_time = []

    with open(thar_list, 'r') as f:
        for line in f:
            name = line.strip()
            try:
                prihdr = fits.getheader(name)
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
                print(f"Can't open file {line}")
    f.close()
    thar_time = np.asarray(thar_time)
    return (thar_name[np.argmin(thar_time)], thar_time[np.argmin(thar_time)])
