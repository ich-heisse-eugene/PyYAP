Pkg_path = "/Volumes/DATA0/work/PyYAP"
# Import general modules
import time
from datetime import datetime, date, time
import os
import glob
import shutil
from sys import path, exit
import argparse
import logging
import multiprocessing as mp
# Below we load pipeline specific modules
from pyyap_logging import listener as log_listener
from lister import lister
from trimmer import trimmer
from fixpix import fixpix
from medianer import medianer
from list_subtractor import list_subtractor
from list_dark_subtractor import dark_subtractor
from thar_combiner import thar_combiner
from tracer import order_tracer, locate_orders
from remap_orders import remap_existing
from sl_remover import sl_subtract
from list_cosmo_cleaner import list_cosmo_cleaner
from extractor import fox, extract_multi
from blz_correct import remove_blz
from thar_reident import thar_auto
from thar_manager import thar_manager
from disp_add import disp_add
from get_sp_resolv import get_sp_resolv
from update_snr import update_snr

## Disable warnings
import warnings
warnings.simplefilter("ignore")
warnings.filterwarnings("ignore")

Pkg_path = os.path.realpath(Pkg_path)
ver = "2025.05"

#####################################################################################
## Parameters
devices = ['mres', 'umres'] # List of valid devices

#####################################################################################
## Let's get started
def S_EX(conf):
    ### Setting up logging facilities
    ##
    m = mp.Manager()
    queue = m.Queue(-1)
    listener = mp.Process(target=log_listener, args=(queue, os.path.join(conf['Path2Data'], 'pyyap_journal.log'), logging.INFO,))
    listener.start()
    logging.getLogger("matplotlib").propagate = False
    #

    start = datetime.now()
    print('Started at:', start.time())

    queue.put((logging.INFO, "PyYAP - PYthon Yet Another Pipeline\n by Eugene Semenko and Vadim Krushinsky"))
    queue.put((logging.INFO, f"Started at: {start.time()}"))
    queue.put((logging.INFO, "=========================================="))
    queue.put((logging.INFO, f"Data in: \'{conf['Path2Data']}\'"))
    queue.put((logging.INFO, f"Device of origin: {conf['device']}"))
    queue.put((logging.INFO, "Current configuration:"))

    print("==========================================")
    print("PyYAP - PYthon Yet Another Pipeline\n by Eugene Semenko and Vadim Krushinsky")
    print("Reduction of echelle spectra")
    print(f"Data in: \'{conf['Path2Data']}\'")
    print(f"Device of origin: {conf['device']}")
    print("Current configuration:")
    for keys, values in conf.items():
        print(f"{keys} : {values}")
        queue.put((logging.INFO, f"{keys} : {values}"))

    Path2Data = conf['Path2Data']
    Path2Raw = os.path.join(Path2Data, 'raw')

    ### Make directory for backing up the raw frames
    ##
    queue.put((logging.INFO, "Make directory for backing up the raw frames"))
    if not os.path.exists(Path2Raw):
        os.makedirs(Path2Raw)
    queue.put((logging.INFO, "...done"))
    #

    ### Make temporary directory for storing raw frames
    ##
    queue.put((logging.INFO, "Create temporary directory for storing raw frames"))
    Path2Temp = os.path.join(Path2Data, 'temp')
    if not os.path.exists(Path2Temp):
        os.makedirs(Path2Temp)
    queue.put((logging.INFO, "...done"))
    #

    ### Check directories and then create four lists with biases, thars, flats, and objects
    ##
    queue.put((logging.INFO, "Now create four lists with biases, thars, flats, and objects"))
    files = lister(Path2Data, Path2Raw, Path2Temp, queue)
    if files != None:
        print ('Lists created')
        queue.put((logging.INFO, " ... lists created"))
    #

    ### Trim overscan and flip spectra
    ##
    queue.put((logging.INFO, "Trim overscan and flip images"))
    flip = conf['flip']
    area = list(map(int, conf['area'].strip('][').split(',')))
    print("Trim overscan and flip image")

    ### Fix CCD cosmetics if necessary
    ##
    if conf['mask'] != 'None':
        queue.put((logging.INFO, "First, fix cosmetics"))
        print("First, fix cosmetics")
        mask = conf['mask']
        status = fixpix(Path2Data, os.path.join(Path2Temp, 'flat_list.txt'), mask, area, flip)
        print(f"= flat frames {status}")
        queue.put((logging.INFO, f"= flat frames {status}"))
        status = fixpix(Path2Data, os.path.join(Path2Temp, 'thar_list.txt'), mask, area, flip)
        print(f"= ThAr frames {status}")
        queue.put((logging.INFO, f"= ThAr frames {status}"))
        status = fixpix(Path2Data, os.path.join(Path2Temp, 'obj_list.txt'), mask, area, flip)
        print(f"= scientific frames {status}")
        queue.put((logging.INFO, f"= scientific frames {status}"))
        queue.put((logging.INFO, "... done"))
    #

    trimmer_data = trimmer(Path2Data, os.path.join(Path2Temp, 'bias_list.txt'), area, flip)
    print(f"= biases {trimmer_data}")
    queue.put((logging.INFO, f" = biases {trimmer_data}"))

    trimmer_data = trimmer(Path2Data, os.path.join(Path2Temp, 'flat_list.txt'), area, flip)
    print(f"= flats {trimmer_data}")
    queue.put((logging.INFO, f"= flats {trimmer_data}"))

    trimmer_data = trimmer(Path2Data, os.path.join(Path2Temp, 'thar_list.txt'), area, flip)
    print(f"= ThAr {trimmer_data}")
    queue.put((logging.INFO, f"= ThAr {trimmer_data}"))

    trimmer_data = trimmer(Path2Data, os.path.join(Path2Temp, 'obj_list.txt'), area, flip)
    print(f"= scientific frames {trimmer_data}")
    queue.put((logging.INFO, f"= scientific frames {trimmer_data}"))
    queue.put((logging.INFO, "... done"))
    #

    ### Remove cosmic hits from images using Pieter G. van Dokkum's L.A.Cosmic Laplacian Cosmic Ray Identification
    ##
    print("Remove cosmic hits from scientific images")
    queue.put((logging.INFO, "Remove cosmic hits from scientific images"))
    status = list_cosmo_cleaner(Path2Temp, 'obj_list.txt', 'obj_CRR_cleaned_list.txt', conf, queue)
    print(f"= objects {status}")
    queue.put((logging.INFO, f"... objects {status}"))
    #

    ### Bias correction
    ##  Create master bias and remove original files
    s_bias_name = conf['s_bias_name']
    print('Create master bias')
    queue.put((logging.INFO, "Create master bias"))
    sbias_data = medianer(Path2Data, os.path.join(Path2Temp, 'bias_list.txt'), s_bias_name)
    print(f"Master bias statistic: Mean = {sbias_data[0]:.2f} Median = {sbias_data[1]:.2f} Sigma = {sbias_data[2]:.2f}")
    queue.put((logging.INFO, f"Master bias statistic: Mean = {sbias_data[0]:.2f} Median = {sbias_data[1]:.2f} Sigma = {sbias_data[2]:.2f}"))
    # Subtract master bias from flats
    print('Subtract master bias from flats')
    queue.put((logging.INFO, "Subtract master bias from flats"))
    status = list_subtractor(os.path.join(Path2Temp, 'flat_list.txt'), os.path.join(Path2Data, s_bias_name), 'Bias', conf)
    print (f"= flats {status}")
    queue.put((logging.INFO, f"Flats: {status}"))
    # Subtract master bias from thars
    print('Subtract master bias from comparison spectra')
    queue.put((logging.INFO, "Subtract master bias from comparison spectra"))
    status = list_subtractor(os.path.join(Path2Temp, 'thar_list.txt'), os.path.join(Path2Data, s_bias_name), 'Bias', conf)
    print(f"= ThAr spectra {status}")
    queue.put((logging.INFO, f"ThAr lamp: {status}"))
    # Subtract master bias from objects
    print('Subtract master bias from objects')
    queue.put((logging.INFO, "Subtract master bias from objects"))
    status = list_subtractor(os.path.join(Path2Temp, 'obj_CRR_cleaned_list.txt'), os.path.join(Path2Data, s_bias_name), 'Bias', conf)
    print(f" = objects {status}")
    queue.put((logging.INFO, f"Objects: {status}"))
    queue.put((logging.INFO, "... done"))
    #

    ### Process dark frames if provided
    ##
    if os.path.isfile(os.path.join(Path2Temp, 'dark_list.txt')):
        print('Start subtracting dark current from a series')
        queue.put((logging.INFO, "Start subtracting dark current from a series"))
        status = dark_subtractor(Path2Data, Path2Temp, ['flat_list.txt', 'thar_list.txt', 'obj_CRR_cleaned_list.txt'], area, flip, s_bias_name, queue)
        print(f"Dark subtraction: {status}")
        queue.put((logging.INFO, f"Dark subtraction: {status}"))
    #

    ### Create master flat and remove old flats
    ##
    s_flat_name = conf['s_flat_name']
    print('Create master flat')
    queue.put((logging.INFO, "Create master flat"))
    sflat_data = medianer(Path2Data, os.path.join(Path2Temp, 'flat_list.txt'), s_flat_name)
    print("... created")
    print(f"Master flat statistics: Mean = {sflat_data[0]:.2f} Median = {sflat_data[1]:.2f} Sigma = {sflat_data[2]:.2f}")
    queue.put((logging.INFO, "... created"))
    queue.put((logging.INFO, f"Master flat statistics: Mean = {sflat_data[0]:.2f} Median = {sflat_data[1]:.2f} Sigma = {sflat_data[2]:.2f}"))
    #

    ### Define some parameters required for order tracing
    #
    slice_half_width = float(conf['slice_half_width'])
    step = int(conf['step'])
    min_height = float(conf['min_height'])
    aperture = float(conf['aperture'])
    view = eval(conf['view'])
    adaptive = eval(conf['adaptive'])
    #

    ### Create master ThAr and delete old files
    ##
    print('Group ThArs by their exposure time and create master ThArs')
    queue.put((logging.INFO, "Group ThArs by their exposure time and create master ThArs"))
    thars = thar_combiner(Path2Data, os.path.join(Path2Temp, 'thar_list.txt'), queue)
    print(f"Master ThArs saved in {thars}")
    queue.put((logging.INFO, f"Master ThArs saved in {thars}"))
    #

    ### Create list of files for averaging to produce the files with distinctive orders
    ##
    print("Orders detection")
    queue.put((logging.INFO, "Orders detection"))
    if 's_ordim_method' in conf:      #  Valid methods: 'hybrid', 'flats', 'objects', and 'remap'
        s_ordim_method = conf['s_ordim_method'].rstrip()
    else:
        s_ordim_method = "hybrid"
    print(f"Method: {s_ordim_method}")
    queue.put((logging.INFO, f"Method: {s_ordim_method}"))
    mapped = False
    if s_ordim_method == "remap":
        mapped = remap_existing(Path2Data, Pkg_path, thars, conf, queue)
    if not mapped:
        mapped = locate_orders(Path2Data, Path2Temp, conf, s_ordim_method, queue)
    if not mapped:
        print("There are no defined orders. Exitting.")
        queue.put((logging.INFO, "There are no defined orders. Exitting."))
        exit(1)
    #

    ### Remove scattered light
    ##
    subtract = eval(conf['subtract'])
    if subtract:
        queue.put((logging.INFO, "Remove scattered light"))
        print("Remove scattered light")
        flat_name, list_name = sl_subtract(Path2Data, Path2Temp, conf, subtract, queue)
    else:
        flat_name = os.path.join(Path2Data, 's_flat.fits')
        list_name = os.path.join(Path2Temp, 'obj_CRR_cleaned_list.txt')
    #

    ### Data extraction
    ##  Extract flat
    print("Data extraction\n= master flat")
    queue.put((logging.INFO, "Data extraction\n= master flat"))
    aperture = float(conf['aperture'])
    ex_type = conf['ex_type']
    ap_file = os.path.join(Path2Temp, 'traces.txt')
    if ex_type !='FOX':
        status, status_err = fox(flat_name, flat_name, ap_file, ex_type, aperture, queue)
        print(f"Extracted s_flat spectrum saved in {status}")
        queue.put((logging.INFO, f"Extracted s_flat spectrum saved in {status}"))
    ##  Extract objects
    print("= objects")
    queue.put((logging.INFO, "= objects"))
    out_list_name = os.path.join(Path2Temp, 'obj_extracted.txt')
    out_err_list_name = os.path.join(Path2Temp, 'err_extracted.txt')
    extract_multi(Path2Data, Path2Temp, list_name, out_list_name, out_err_list_name, conf, flat_name, queue)
    ## Extract ThAr
    print("= ThAr")
    queue.put((logging.INFO, "= ThAr"))
    list_name = os.path.join(Path2Temp, 's_thar_list.txt')
    out_list_name = os.path.join(Path2Temp, 's_thar_extracted.txt')
    f_out=open(out_list_name, 'a')
    with open(list_name, 'r') as f:
        for line in f:
            name = line.strip()
            print(f"flat_name = {flat_name}, ap_file = {ap_file}, aperture = {aperture}")
            queue.put((logging.INFO, f"flat_name = {flat_name}, ap_file = {ap_file}, aperture = {aperture}"))
            status, _ = fox(name, flat_name, ap_file, 'APEX', aperture, queue) # Note the fixed method of extraction
            print(f"Extracted spectrum saved in {status}")
            queue.put((logging.INFO, f"Extracted spectrum saved in {status}"))
            shutil.move(os.fspath(name), os.fspath(Path2Temp))
            print(status, file=f_out)
    f.close()
    f_out.close()
    os.remove(list_name)
    #

    ### Stop processing here in case of extraction without wavelength calibration
    ##
    if not eval(conf['calibrate']):
        try:
            os.mkdir(os.path.join(Path2Data, 'Reduced'))
        except OSError as e:
            print(f"Exception: {e}")
        finally:
            source = os.listdir(Path2Data)
            for files in source:
                if files.endswith('_ec.fits'):
                    shutil.move(os.fspath(files), os.fspath(os.path.join(Path2Data, 'Reduced')))
                if files.endswith('_err.fits'):
                    shutil.move(os.fspath(files), os.fspath(os.path.join(Path2Data, 'Reduced')))
                if files.endswith('.pdf'):
                    shutil.move(os.fspath(files), os.fspath(os.path.join(Path2Data, 'Reduced')))
                if files.endswith('.log'):
                    shutil.move(os.fspath(files), os.fspath(os.path.join(Path2Data, 'Reduced')))

        if eval(conf['strip']):
            print("Cleaning catalogues...")
            queue.put((logging.INFO, "Cleaning catalogues..."))
            shutil.rmtree(os.path.join(Path2Data, 'temp'),  ignore_errors=True)
            for filepath in glob.iglob(os.path.join(Path2Data, '*.*')):
                os.remove(filepath)

        end = datetime.now()
        print(f"Finished at: {end.time()}")
        queue.put((logging.INFO, f"Finished at: {end.time()}"))
        i, d = divmod((end-start).seconds/60, 1)
        print(f"Duration (m:s): {i:3.0f}:{int(d*60)}")
        queue.put((logging.INFO, f"Duration (m:s): {i:3.0f}:{int(d*60)}"))
        return None
    #

    ### Otherwise search and identify lines in ThAr
    ##
    list_name = os.path.join(Path2Temp, 's_thar_extracted.txt')
    # copy instrument-related files with existing calibrations to the working directory
    shutil.copy2(os.path.join(Pkg_path, 'devices', conf['device'], conf['device']+'_thar.dat'), os.path.join(Path2Data, 'thar.dat'))
    shutil.copy2(os.path.join(Pkg_path, 'devices', conf['device'], conf['device']+'_thar_last.dat'), os.path.join(Path2Data, 'thar_last.dat'))
    shutil.copy2(os.path.join(Pkg_path, 'devices', conf['device'], conf['device']+'_thar_last.fits'), os.path.join(Path2Data, 'thar_last.fits'))
    if not os.path.exists(os.path.join(Path2Data, 'Old_Thars')):
        os.makedirs(os.path.join(Path2Data, 'Old_Thars'))

    with open(list_name, 'r') as f:
        anr = float(conf['anr'])
        xord = int(conf['xord'])
        yord = int(conf['yord'])
        for line in f:
            name = line.strip()
            print(f"Search WL solution for: {name}")
            queue.put((logging.INFO, f"Search WL solution for: {name}"))
            thar_auto(Path2Data, name, anr, xord, yord, view, queue)
            ## Save the new WL solution to archive
            dt=datetime.now()
            dt=dt.strftime("%y-%m-%dT%H-%M")
            shutil.copy2(name, os.path.join(Path2Data, 'Old_Thars', dt + '_thar.fits'))
            shutil.copy2(name.replace('_ec', '_err'), os.path.join(Path2Data, 'Old_Thars', dt + '_thar_err.fits'))
            shutil.copy2(os.path.splitext(name)[0] + "_disp.txt", os.path.join(Path2Data, 'Old_Thars', dt + '_thar_disp.txt'))
            shutil.copy2(os.path.splitext(name)[0] + "_features.txt", os.path.join(Path2Data, 'Old_Thars', dt + '_thar_features.txt'))
    f.close()
    #

    shutil.move(os.fspath(flat_name), os.fspath(Path2Temp))
    os.remove(os.path.join(Path2Data, 'thar_last.dat'))
    os.remove(os.path.join(Path2Data, 'thar.dat'))
    os.remove(os.path.join(Path2Data, 'thar_last.fits'))
    dt=datetime.now()
    dt=dt.strftime("%y-%m-%dT%H-%M")
    shutil.copy2(os.path.join(Path2Temp, 'traces.txt'), os.path.join(Path2Data, 'Old_Thars', dt + '_traces.txt'))

    ### Search for a calibration nearest to a moment of exposure and then apply the WL solution
    ##
    obj_list = os.path.join(Path2Temp, 'obj_extracted.txt')
    obj_err_list = os.path.join(Path2Temp, 'err_extracted.txt')
    thar_list = os.path.join(Path2Temp, 's_thar_extracted.txt')
    # ThAr
    with open(thar_list, 'r') as f:
        for line in f:
            name = line.strip()
            solution = thar_manager(name, thar_list, queue)
            print(f"{name}: nearest thar is {solution[0].split(os.sep)[-1]}, difference: {solution[1]} sec")
            queue.put((logging.INFO, f"{name}: nearest thar is {solution[0].split(os.sep)[-1]}, difference: {solution[1]} sec"))
            name_cal = disp_add(name, solution[0], view, queue)
            get_sp_resolv(name_cal, queue)
    f.close()
    # Variances
    with open(obj_err_list, 'r') as f:
        for line in f:
            name = line.strip()
            print(name)
            solution = thar_manager(name, thar_list, queue)
            print(f"Nearest thar is {solution[0].split(os.sep)[-1]}, difference: {solution[1]} sec")
            queue.put((logging.INFO, f"Nearest thar is {solution[0].split(os.sep)[-1]}, difference: {solution[1]} sec"))
            name_cal = disp_add(name, solution[0], view, queue)
    f.close()
    # Objects
    with open(obj_list, 'r') as f:
        for line in f:
            name = line.strip()
            print(name)
            solution = thar_manager(name, thar_list, queue)
            print(f"{name}: nearest thar is {solution[0].split(os.sep)[-1]}, difference: {solution[1]} sec")
            queue.put((logging.INFO, f"{name}: nearest thar is {solution[0].split(os.sep)[-1]}, difference: {solution[1]} sec"))
            name_cal = disp_add(name, solution[0], view, queue)
            update_snr(name_cal, queue)
    f.close()
    #

    ### Divide images by flat
    ##
    if ex_type !='FOX':
        name = str(flat_name).replace('.fits', '_ec.fits')
        solution = thar_manager(name, thar_list, queue)
        print(f"{name}: nearest thar is {str(solution[0].split(os.sep)[-1])}, difference: {str(solution[1])} sec")
        queue.put((logging.INFO, f"{name}: nearest thar is {str(solution[0].split(os.sep)[-1])}, difference: {str(solution[1])} sec"))
        name_cal = disp_add(name, solution[0], view, queue)
        print("Flat field normalisation of data ")
        queue.put((logging.INFO, "Flat field normalisation of data "))
        obj_list = os.path.join(Path2Temp, 'obj_extracted.txt')
        s_blaze_name = name_cal
        with open(obj_list, 'r') as f:
            for line in f:
                name = line.strip()
                in_name = str(name).replace('_ec.fits', '_ec_WCS.fits')
                out_name = str(name).replace('_ec.fits', '_blz_ec_WCS.fits')
                status = remove_blz(in_name, s_blaze_name, out_name)
                print(f"Corrected spectrum was saved in {status}")
                queue.put((logging.INFO, f"Corrected spectrum was saved in {status}"))
    #

    ### Finishing up
    ##
    try:
        if not os.path.exists(os.path.join(Path2Data, 'Reduced')):
            os.mkdir(os.path.join(Path2Data, 'Reduced'))
    except OSError as e:
        print(f"Exception: {e}")
    finally:
        source = os.listdir(Path2Data)
        for files in source:
            if os.path.isfile(files) and not files.endswith('.log'):
                if files.endswith('_WCS.fits'):
                    shutil.copy2(os.fspath(files), os.fspath(os.path.join(Path2Data, 'Reduced')))
                elif files.endswith('.pdf'):
                    shutil.copy2(os.fspath(files), os.fspath(os.path.join(Path2Data, 'Reduced')))
                elif files.endswith('.log'):
                    shutil.copy2(os.fspath(files), os.fspath(os.path.join(Path2Data, 'Reduced')))
                else:
                    shutil.copy2(os.fspath(files), os.fspath(os.path.join(Path2Data, 'temp')))
                os.remove(files)

    if eval(conf['strip']):
        print("Cleaning up catalogues...")
        queue.put((logging.INFO, "Cleaning up catalogues"))
        shutil.rmtree(os.path.join(Path2Data, 'temp'),  ignore_errors=True)
        shutil.rmtree(os.path.join(Path2Data, 'Old_Thars'),  ignore_errors=True)

    end = datetime.now()
    print(f"Finished at: {end.time()}")
    queue.put((logging.INFO, f"Finished at: {end.time()}"))
    i, d = divmod((end-start).seconds/60, 1)
    print(f"Duration (m:s): {i:3.0f}:{int(d*60)}")
    queue.put((logging.INFO, f"Duration (m:s): {i:3.0f}:{int(d*60)}"))
    queue.put_nowait(None)
    listener.join()
    shutil.copy2(os.fspath(os.path.join(Path2Data, 'pyyap_journal.log')), os.fspath(os.path.join(Path2Data, 'Reduced')))
    os.remove(os.fspath(os.path.join(Path2Data, 'pyyap_journal.log')))
    return None

## Main programme  ##
if __name__ == "__main__":
    conf = {}
    parser = argparse.ArgumentParser()
    parser.add_argument("cwd", help="Current work directory")
    parser.add_argument("--device", help="Specify the origin of input data", type=str, default="mres")
    parser.add_argument("--extractonly", help="Extraction without calibration", action="store_true")
    parser.add_argument("--sl", help="Turn on/off the subtraction of scattered light [True/False]", default=None)
    parser.add_argument("--method", help="Method of extraction [FOX, APEX, PSFEX]", default=None)
    parser.add_argument("--strip", help="Clean supplement materials [True/False]", default=None)
    parser.add_argument("--view", help="Turn on/of the visualization of the results [True/False]", default=None)
    args = parser.parse_args()

    Data_path = os.path.realpath(args.cwd.strip())
    conf['Path2Data'] = Data_path
    conf['device'] = args.device.strip()
    if os.path.isfile(os.path.join(Data_path, 'names.txt')):
        from check_header import fill_headers
        print(f"Correcting FITS headers in {Data_path} \r")
        fill_headers(os.path.join(Data_path, 'names.txt'), conf['device'], ver)
    if not os.path.isfile(os.path.join(Data_path, args.device + '.conf')):
        shutil.copy(os.path.join(Pkg_path, 'devices', args.device, args.device+'.conf'), Data_path)
    with open(os.path.join(Data_path, args.device + '.conf')) as f:
        for line in f.readlines():
            if not line.strip().startswith("#"):
                (key, val) = line.split('=')
                conf[key.strip()] = val.strip()

    if args.sl != None and args.sl in ['True', 'False']:
        conf['subtract'] = args.sl
    if args.method != None and args.method in ['FOX', 'APEX', 'PSFEX']:
        conf['ex_type'] = args.method
    if args.strip != None and args.strip in ['True', 'False']:
        conf['strip'] = args.strip
    if args.view != None and args.view in ['True', 'False']:
        conf['view'] = args.view
    if args.extractonly:
        conf['calibrate'] = 'False'
    if conf['mask'].strip() != 'None':
        if not os.path.isfile(os.path.join(Data_path, conf['mask'].strip())):
            conf['mask'] = os.path.join(Pkg_path, 'devices', args.device, conf['mask'].strip())
    ## Main entrance point
    S_EX(conf)
    exit(0)
