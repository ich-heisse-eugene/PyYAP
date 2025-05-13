from astropy.io import fits
import os
import shutil
import logging

##################################################################
def lister(data_dir, raw_dir, temp_dir, queue):
    file_list = []
    darks = False
    dir_content = os.listdir(data_dir)

    #create list of fits-files
    for ii in range (0,len(dir_content)):
        if dir_content[ii].count('.fit')==1 or dir_content[ii].count('.fts')==1 or dir_content[ii].count('.fits')==1:
            file_list.append(os.path.join(data_dir, dir_content[ii]))

    counter = len(file_list)
    if len(file_list) > 0:
        bf = open(os.path.join(temp_dir, 'bias_list.txt'), 'w')
        df = open(os.path.join(temp_dir, 'dark_list.txt'), 'w')
        ff = open(os.path.join(temp_dir, 'flat_list.txt'), 'w')
        tf = open(os.path.join(temp_dir, 'thar_list.txt'), 'w')
        of = open(os.path.join(temp_dir, 'obj_list.txt'), 'w')
        #create lists
        for ii in range (0,len(file_list)):
            try:
                shutil.copy2(file_list[ii], raw_dir)
                hdulist = fits.open(file_list[ii], mode='update')
                prihdr = hdulist[0].header
                if 'IMAGETYP' not in prihdr:
                    if os.fspath(file_list[ii]).lower().find('bias') != -1:
                        im_type = 'BIAS'
                    if os.fspath(file_list[ii]).lower().find('dark') != -1:
                        im_type = 'DARK'
                    elif os.fspath(file_list[ii]).lower().find('flat') != -1:
                        im_type = 'FLAT'
                    elif os.fspath(file_list[ii]).lower().find('thar') != -1:
                        im_type = 'THAR'
                    elif os.fspath(file_list[ii]).lower().find('obj') != -1:
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
                elif im_type.lower().find('dark') != -1:
                    if 'IMAGETYP' in prihdr:
                        prihdr['IMAGETYP'] = 'DARK'
                    else:
                        prihdr.set('IMAGETYP', 'DARK', 'Dark current frame')
                    if darks == False: darks = True
                    print(file_list[ii], file=df)
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
                queue.put((logging.INFO, "Err: Can't open file: {file_list[ii]}"))
                pass
        bf.close()
        df.close()
        ff.close()
        tf.close()
        of.close()
        if not darks:
            os.remove(os.path.join(temp_dir, 'dark_list.txt'))
        if counter==0:
            queue.put((logging.INFO, "OK"))
            return ('OK')
        else:
            print('Something is wrong with files')
            queue.put((logging.INFO, "Err: Something is wrong with files"))
            return(None)
    else:
        print('Directory is empty')
        queue.put((logging.INFO, "Err: Directory is empty"))
        return(None)
