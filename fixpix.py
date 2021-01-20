import astropy.io.fits as pyfits
import numpy as np
from scipy.interpolate import splev, splrep
import os
import copy
import shutil

"""
Created by Eugene Semenko 29.07.2020
This function fixes cosmetics of CCD images according to the mask
provided by user.
"""

##################################################################
def replace_bad_data(img, idx, idy, columns):
    npix = 5
    if columns:
        for y in idy:
            dy = np.arange(y[0]-npix, y[-1]+npix+1, 1, dtype=int)
            dy_interp = np.arange(y[0], y[-1]+1, 1, dtype=int)
            img_idx = np.vstack((np.arange(y[0]-npix, y[0], 1, dtype=int), np.arange(y[-1]+1, y[-1]+npix+1, 1, dtype=int))).flatten()
            for x in idx:
                for coord in range(x[0], x[-1]+1, 1):
                    spl = splrep(img_idx, img[coord, img_idx])
                    img[coord, dy_interp] = splev(dy_interp, spl)
    else:
        for x in idx:
            dx = np.arange(x[0]-npix, x[-1]+npix+1, 1, dtype=int)
            dx_interp = np.arange(x[0], x[-1]+1, 1, dtype=int)
            img_idx = np.vstack((np.arange(x[0]-npix, x[0], 1, dtype=int), np.arange(x[-1]+1, x[-1]+npix+1, 1, dtype=int))).flatten()
            for y in idy:
                for coord in range(y[0], y[-1]+1, 1):
                    spl = splrep(img_idx, img[img_idx, coord])
                    img[dy_interp, coord] = splev(dy_interp, spl)
    return img


def fixpix(dir_name, list_name, mask_file, area, flip):
    columns = True
    # Define the map of bad pixels
    try:
        mhdu = pyfits.open(mask_file)
        mask = mhdu[0].data[0].copy()
        mhdu.close()
    finally:
        pass
        # print(f"Mask file is {mask_file}")

    # trim and flip mask
    if mask.shape[1]>(area[1]-area[0]) and mask.shape[0]>(area[3]-area[2]):
        trimmed_mask = copy.copy(mask[area[2]:area[3],area[0]:area[1]])
        if flip == 'X':
            flipped_mask = np.flip(trimmed_mask, 1)
        elif flip == 'Y':
            flipped_mask = np.flip(trimmed_mask, 0)
        elif flip == 'XY':
            flipped_mask = np.flip(trimmed_mask)
        else:
            flipped_mask = trimmed_mask.copy()
    mask = flipped_mask

    midx = np.where(mask == 0)
    nmx = np.unique(midx[0])
    nmy = np.unique(midx[1])
    if len(nmy) < len(nmx):
        idy = np.split(nmy, np.where(np.diff(nmy) != 1)[0] + 1)
        idx = []
        for col in idy:
            indices = np.where(mask[:,col] == 0)
            idx.append(np.array([np.min(indices), np.max(indices)], dtype=int))
    else:
        columns = False
        idx = np.split(nmx, np.where(np.diff(nmx) != 1)[0] + 1)
        idy = []
        for row in idx:
            indices = np.where(mask[row,:] == 0)
            idy.append(np.array([np.min(indices), np.max(indices)], dtype=int))

    with open(dir_name.joinpath(list_name), 'r') as f:
        for line in f:
            name = line.strip()
            try:
                hdulist = pyfits.open(name, mode = 'update')
                data = hdulist[0].data.copy()
                prihdr = hdulist[0].header
                hdulist.close()

                img_fixed = replace_bad_data(data, idx, idy, columns)
                hdulist[0].data = img_fixed
                prihdr['HISTORY'] = 'CCD cosmetics fixed'
                try:
                    hdulist.writeto(name, overwrite=True)
                except IOError:
                    print ("ERROR: Can't write file '%s'" % name)
            except IOError:
                print ("Can't open file:", name)
    f.close()

    return ("Fixed")
