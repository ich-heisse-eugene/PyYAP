from LA_Cosmic import detCos
import shutil
import os
import logging
import multiprocessing as mp

def process(name, queue):
    out_name = os.path.splitext(name)[0] + '_CRR.fits'
    detCos(image=name,  out_clean=out_name, verbose=True)
    print(f"{out_name} created")
    queue.put((logging.INFO, f"{out_name} created"))
    return out_name

def errclbk(error):
    print(error)
    return None

def list_cosmo_cleaner(dir_name, list_name, out_list_name, conf, queue):
    print("Starting")
    queue.put((logging.INFO, "Starting"))
    list_out = []
    with open(os.path.join(dir_name, list_name), 'r') as f:
        nCPUs = os.cpu_count()
        if 'threading' in conf and eval(conf['threading']) and nCPUs > 2:
            queue.put((logging.INFO, "Multiprocessing is ON"))
            proc_args = [(line.strip(), queue) for line in f]
            with mp.Pool(processes=nCPUs) as pool:
                res_async = pool.starmap_async(process, proc_args, chunksize=nCPUs, error_callback=errclbk)
                list_out.extend(res_async.get())
        else:
            for line in f:
                res_mono = process(line.strip(), queue)
                list_out.extend([res_mono])

    with open(os.path.join(dir_name, out_list_name), 'a') as f_out:
        for item in list_out:
            print(item, file=f_out)

    if os.path.isfile(os.path.join(dir_name, out_list_name)):
        print(f"File names saved in {os.path.join(dir_name, out_list_name)}")
        queue.put((logging.INFO, f"File names saved in {os.path.join(dir_name, out_list_name)}"))
    else:
        print("Cleaned files haven't been saved")
    return "cleaned"
