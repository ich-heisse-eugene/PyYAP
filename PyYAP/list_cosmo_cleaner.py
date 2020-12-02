from LA_Cosmic import detCos
import shutil
import os
import logging

##################################################################

def list_cosmo_cleaner(dir_name, list_name, out_list_name):

    f_out=open(dir_name+'/'+out_list_name, 'a')

    with open(dir_name+'/'+list_name, 'r') as f:
        for line in f:
            name = line.strip()
            out_name = os.path.splitext(name)[0] + '_CRR.fits'
            detCos(image=name,  out_clean=out_name)
            f_out.write(out_name + '\n')
            print(out_name)
            logging.info(f"{out_name}")
            # shutil.move(name, dir_name +'/temp')
            print()

    f.close()
    f_out.close()

    # os.remove(dir_name+list_name)

    print('File names saved in ', dir_name+out_list_name)
    return ("Cleaned")

#################################test######################################
# dir_name = '/home/eugene/pipelines/PyYAP/20200112/temp'
# list_name = '/obj_list.txt'
# out_list_name = '/obj_CRR_cleaned_list.txt'
# ##
# list_cosmo_cleaner(dir_name, list_name, out_list_name)
