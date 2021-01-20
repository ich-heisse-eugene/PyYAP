import astropy.io.fits as pyfits
import os
import shutil
import logging

##################################################################
def lister(data_dir, raw_dir, temp_dir):
    file_list = []
    dir_content = os.listdir(data_dir)

    #create list of fits-files
    for ii in range (0,len(dir_content)):
        if dir_content[ii].count('.fit')==1 or dir_content[ii].count('.fts')==1 or dir_content[ii].count('.fits')==1:
            file_list.append(data_dir.joinpath(dir_content[ii]))

    counter = len(file_list)
    if len(file_list) > 0:
        bf = open(temp_dir.joinpath('bias_list.txt'), 'w')
        ff = open(temp_dir.joinpath('flat_list.txt'), 'w')
        tf = open(temp_dir.joinpath('thar_list.txt'), 'w')
        of = open(temp_dir.joinpath('obj_list.txt'), 'w')
        #create lists
        for ii in range (0,len(file_list)):
            try:
                shutil.copy2(file_list[ii], raw_dir)
                hdulist = pyfits.open(file_list[ii], mode='update')
                prihdr = hdulist[0].header
                if 'IMAGETYP' not in prihdr:
                    if file_list[ii].lower().find('bias') != -1:
                        im_type = 'BIAS'
                    elif file_list[ii].lower().find('flat') != -1:
                        im_type = 'FLAT'
                    elif file_list[ii].lower().find('thar') != -1:
                        im_type = 'THAR'
                    elif file_list[ii].lower().find('obj') != -1:
                        im_type = 'OBJ'
                else:
                    im_type = prihdr['IMAGETYP']
                if im_type.lower().find('bias') != -1:
                    print(file_list[ii], file=bf)
                    counter = counter-1
                elif im_type.lower().find('flat') != -1:
                    if 'IMAGETYP' in prihdr:
                        prihdr['IMAGETYP'] = 'FLAT'
                    else:
                        prihdr.set('IMAGETYP', 'FLAT', 'Flat field spectrum')
                    print(file_list[ii], file=ff)
                    counter = counter-1
                elif im_type.lower().find('thar') != -1:
                    if 'IMAGETYP' in prihdr:
                        prihdr['IMAGETYP'] = 'THAR'
                    else:
                        prihdr.set('IMAGETYP', 'THAR', 'ARC spectrum')
                    print(file_list[ii], file=tf)
                    counter = counter-1
                elif im_type.lower().find('obj') != -1:
                    if 'IMAGETYP' in prihdr:
                        prihdr['IMAGETYP'] = 'OBJ'
                    else:
                        prihdr.set('IMAGETYP', 'OBJ', 'Scientific spectrum')
                    print(file_list[ii], file=of)
                    counter = counter-1
                hdulist[0].header = prihdr
                hdulist.flush()
                hdulist.close()
            except IOError:
                print("Can't open file: ", file_list[ii])
                logging.error("Can't open file: {file_list[ii]}")
                pass
        bf.close()
        ff.close()
        tf.close()
        of.close()
        if counter==0:
            logging.info('OK')
            return ('OK')
        else:
            print('Some files are wrong')
            logging.warning('Some files are wrong')
            return(None)
    else:
        print('Directory is empty')
        return(None)
