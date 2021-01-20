from LA_Cosmic import detCos
import shutil
import os
import logging

##################################################################

def list_cosmo_cleaner(dir_name, list_name, out_list_name):

    f_out=open(dir_name.joinpath(out_list_name), 'a')

    with open(dir_name.joinpath(list_name), 'r') as f:
        for line in f:
            name = line.strip()
            out_name = os.path.splitext(name)[0] + '_CRR.fits'
            detCos(image=name,  out_clean=out_name)
            print(out_name, file=f_out)
            print(out_name)
            logging.info(f"{out_name}")
            print()

    f.close()
    f_out.close()

    print('File names saved in ', dir_name.joinpath(out_list_name))
    return("Cleaned")
