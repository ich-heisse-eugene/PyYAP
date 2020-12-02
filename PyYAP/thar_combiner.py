import time
import astropy.io.fits as pyfits
import numpy as np
from medianer import medianer
import os
import logging
##import platform

def thar_combiner(dir_name, thar_list):
    thar_name = []
    thar_time = []

    with open(thar_list, 'r') as f:
        for line in f:
            name = line.strip()
            prihdr = pyfits.getheader(name)
            thar_name.append(name.split('/')[-1])
            print(prihdr['DATE'])
            thar_time.append(time.mktime(time.strptime(prihdr['DATE'], "%Y-%m-%dT%H:%M:%S")))#.%f
    f.close()
    print()

    _set = np.zeros(len(thar_time))
    thar_all = []

    for ii in range(0, len(thar_time)):
        thar_ansamble = []
        for jj in range (0, len(thar_time)):
            if _set[jj]!=1:
                if abs(thar_time[jj]-thar_time[ii])<600:
                    _set[jj]=1
                    thar_ansamble.append(thar_name[jj])

        if len(thar_ansamble)>0:
            thar_all.append(thar_ansamble)

    tf=open(dir_name+'/temp/s_thar_list.txt', 'w')
    for ii in range(0, len(thar_all)):
        temp_list = open(dir_name+'/temp.txt', 'w')
        local = thar_all[ii]
        for jj in range(0, len(local)):
            temp_list.write(dir_name+'/'+local[jj] + '\n')
            print (local[jj])
        temp_list.close()
        medianer (dir_name, dir_name+'/temp.txt', 's_thar_'+str(ii)+'.fits')
        tf.write(dir_name+'/s_thar_'+str(ii)+'.fits' + '\n')
        print(f"s_thar_{str(ii)}.fits created")
        logging.info(f"s_thar_{str(ii)}.fits created")
        print()
    tf.close()

    os.remove(thar_list)
    os.remove(dir_name+'/temp.txt')

    return(dir_name+'/temp/s_thar_list.txt')
