import astropy.io.fits as pyfits
import os
import numpy
import copy
import shutil

##################################################################
def trimmer(dir_name, list_name, area, flip):
    with open(dir_name.joinpath(list_name), 'r') as f:
        for line in f:
            name = line.strip()
            try:
                hdulist = pyfits.open(name, mode = 'update', do_not_scale_image_data=True)
                prihdr = hdulist[0].header
                data = hdulist[0].data.copy()
                hdulist.close()
                if data.shape[1]>(int(area[1])-int(area[0])) and data.shape[0]>(int(area[3])-int(area[2])):
                    trimmed_data = copy.copy(data[int(area[2]):int(area[3]),int(area[0]):int(area[1])])
                    if flip == 'X':
                        flipped_data = numpy.flip(trimmed_data, 1)
                    elif flip == 'Y':
                        flipped_data = numpy.flip(trimmed_data, 0)
                    elif flip == 'XY':
                        flipped_data = numpy.flip(trimmed_data)
                    else:
                        flipped_data = trimmed_data.copy()
                    hdulist[0].data = flipped_data
                    prihdr['NAXIS1'] = flipped_data.shape[1]
                    prihdr['NAXIS2'] = flipped_data.shape[0]
                    prihdr['HISTORY'] = 'overscan trimmed'
                    prihdr['HISTORY'] = 'Data flipped along '+flip
                    try:
                        hdulist.writeto(name, overwrite=True)
                    except IOError:
                        print(f"ERROR: Can't write file {name}")
                else:
                    print (f"Frame {name} has wrong size")

            except IOError:
                print (f"Can't open file: {name}")
    f.close()
    return("Trimmed")
