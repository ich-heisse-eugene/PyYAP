from astropy.io import fits
import numpy as np
import logging
import os
import shutil

from trimmer import trimmer
from medianer import medianer
from list_subtractor import list_subtractor
from list_cosmo_cleaner import list_cosmo_cleaner

def parse_exposures(Path2Temp, list_names):
    with open(os.path.join(Path2Temp, list_names), 'r') as fo:
        obj_exp = np.array([], dtype=int)
        obj_fn = np.array([])
        for line in fo:
            line = line.rstrip()
            with fits.open(line) as hdu:
                obj_fn = np.append(obj_fn, line)
                hdr = hdu[0].header
                if 'EXPTIME' in hdr:
                    obj_exp = np.append(obj_exp, int(hdr['EXPTIME']))
    return obj_exp, obj_fn

def process_block(Path2Data, Path2Temp, list_names, dark_time_exp, conf, queue):
    obj_exp, obj_fn = parse_exposures(Path2Temp, list_names)
    obj_time_exp = np.unique(obj_exp)
    for t in obj_time_exp:
        if t in dark_time_exp:
            idx = np.where(obj_exp == t)
            np.savetxt(os.path.join(Path2Temp, 's_temp_dark_'+str(t)+'s.txt'), obj_fn[idx], fmt='%s')
            status = list_subtractor(os.path.join(Path2Temp, 's_temp_dark_'+str(t)+'s.txt'), os.path.join(Path2Data, 's_dark_'+str(t)+'s.fits'), 'Dark', conf)
            print(f"Texp = {t} s dark subtraction: {status}")
            queue.put((logging.INFO, f"Texp = {t} s dark subtraction: {status}"))
            os.remove(os.path.join(Path2Temp, 's_temp_dark_'+str(t)+'s.txt'))
            print()
    return None

def dark_subtractor(Path2Data, Path2Temp, lists_names, area, flip, s_bias_name, conf, queue):
    trimmer_data = trimmer(Path2Data, os.path.join(Path2Temp, 'dark_list.txt'), area, flip)
    print("Darks trimmed")
    queue.put((logging.INFO, "Darks trimmed"))
    status = list_cosmo_cleaner(Path2Temp, 'dark_list.txt', 'dark_CRR_cleaned_list.txt', conf, queue)
    print(f"Darks {status}")
    queue.put((logging.INFO, f"Darks {status}"))
    print()
    status = list_subtractor(os.path.join(Path2Temp, 'dark_CRR_cleaned_list.txt'), os.path.join(Path2Data, s_bias_name), 'Bias', conf)
    print(f"Darks: {status}")
    queue.put((logging.INFO, f"Darks: {status}"))
    print()
    dark_exp = np.array([], dtype=int)
    dark_fn = np.array([])
    dark_exp, dark_fn = parse_exposures(Path2Temp, os.path.join(Path2Temp, 'dark_CRR_cleaned_list.txt'))
    dark_time_exp = np.unique(dark_exp)
    for t in dark_time_exp:
        idx = np.where(dark_exp == t)
        if len(idx[0]) >= 2:
            np.savetxt(os.path.join(Path2Temp, 's_dark_'+str(t)+'s.txt'), dark_fn[idx], fmt='%s')
            sdark_data = medianer(Path2Data, os.path.join(Path2Temp, 's_dark_'+str(t)+'s.txt'), os.path.join(Path2Data, 's_dark_'+str(t)+'s.fits'))
            print(f"Statistics for super dark file with Texp={t} s: Mean = {sdark_data[0]:.2f} Median = {sdark_data[1]:.2f} Sigma = {sdark_data[2]:.2f}")
            queue.put((logging.INFO, f"Statistics for super dark file with Texp={t} s: Mean = {sdark_data[0]:.2f} Median = {sdark_data[1]:.2f} Sigma = {sdark_data[2]:.2f}"))
        elif len(idx[0]) == 1:
            np.savetxt(os.path.join(Path2Temp, 's_dark_'+str(t)+'s.txt'), dark_fn[idx], fmt='%s')
            shutil.copy(dark_fn[idx][0], os.path.join(Path2Data, 's_dark_'+str(t)+'s.fits'))

    for item in lists_names:
        process_block(Path2Data, Path2Temp, item, dark_time_exp, conf, queue)
    return 'completed'
