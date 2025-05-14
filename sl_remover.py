from astropy.io import fits
import os
import shutil
import numpy as np
import multiprocessing as mp
import logging

from surf_fit import surf_fit
from scipy.ndimage import gaussian_filter

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from  matplotlib.colors import Normalize as Normalize

import warnings
warnings.simplefilter("ignore")

##################################################################

def read_traces(x_coo, ap_file):
    Y = []
    FWHM = []
    with open(ap_file, 'r') as fp:
        f = fp.readlines()
        p = f[3].strip().rsplit()
        if len(p) == 4:
            n_orders = int(p[0]); poly_order = int(p[1])
            adaptive = eval(p[2]); ap_size = float(p[3])
        else:
            print("Problem of reading the file with traces")
        if adaptive:
            for i in range(4, 4+n_orders):
                p = f[i].strip().rsplit()
                poly_trace_coef = np.asarray(p[3:3+poly_order+1], dtype=float)
                poly_width_coef = np.asarray(p[3+poly_order+1:], dtype=float)
                Y.append(np.polyval(poly_trace_coef, x_coo))
                FWHM.append(np.polyval(poly_width_coef, x_coo))
        else:
            for i in range(4, 4+n_orders):
                p = f[i].strip().rsplit()
                poly_trace_coef = np.asarray(p[3:3+poly_order+1], dtype=float)
                Y.append(np.polyval(poly_trace_coef, x_coo))
                FWHM.append(np.repeat(float(p[3+poly_order+1]), len(x_coo)))
    return np.asarray(Y), np.asarray(FWHM)

def sl_remover(dir_name, Path2Temp, file_name, ap_file, step, x_order, y_order, subtract, view, queue):
    print(f"Scattered light removal for {file_name} started")
    queue.put((logging.INFO, f"Scattered light removal for {file_name} started"))

    #read image file
    with fits.open(os.path.join(dir_name, file_name)) as hdulist:
        arr = hdulist[0].data.copy()
        prihdr = hdulist[0].header


    orders, _ = read_traces(np.arange(arr.shape[1]), ap_file)

    X = []
    #calculate middle line for 2 close orders
    for jj in range(0, orders.shape[0]-1):
        X.append((orders[jj, :]+orders[jj+1, :])/2)
    X = np.asarray(X)   #array with gaps

    background_X = np.arange(step/2, X.shape[1], step) #array with Xcoo of background probe[number of probe along X]

    Y = []
    for ii in range(0, background_X.shape[0]):
        Y.append(X[:, int(background_X[ii])])

    Y = np.asarray(Y).T #array with Ycoo of background probe [number of probe along X, number of probe along Y]
    x=background_X.tolist()
    XX = np.array([x for i in Y[:,0]])

    ##### Debug. Plot the image with points where the background was taken
    # plt.imshow(arr, cmap='gray')
    # # plt.clim(0, 1000)
    # plt.plot(XX, Y, 'r+', ms=0.5)
    # plt.show()

    #get probe of scatter light
    Y[np.where(Y<1)] = 1
    Z = np.zeros_like(Y)
    for ii in range(0, Y.shape[0]):
        for jj in range(0, Y.shape[1]):
            Z[ii, jj] = np.min(arr[int(Y[ii,jj]-2):int(Y[ii,jj]+3), int(background_X[jj])])

    Z = gaussian_filter(Z, sigma=(2.5, 6.5), mode='nearest')
    # # Debug. Surface plot of scattered light
    if view:
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.plot_surface(XX, Y, Z, cmap='plasma', edgecolor='none')
        ax.set_title('Surface plot')
        plt.show()

    # sc_light = median_filter(surf_fit(XX, Y, Z, x_order, y_order, arr.shape[0], arr.shape[1]), size=15)
    sc_light = surf_fit(XX, Y, Z, x_order, y_order, arr.shape[0], arr.shape[1])

    ##test
    if subtract:
        arr_new=arr-sc_light
    else:
        arr_new=arr
    Y = np.transpose(Y)
    Z_new = np.zeros_like(Y)
    for ii in range(0, Y.shape[0]):
        for jj in range(0, Y.shape[1]):
            Z_new[ii, jj]=np.median(arr_new[int(Y[ii,jj]-1) : int(Y[ii,jj]+1), int(background_X[ii]-step/2):int(background_X[ii]+step/2)])
    Z_new = np.transpose(Z_new)

    print(f"File {file_name} --- {len(orders)} orders\nUse {Y.shape[1]} x {Y.shape[0]} points for scattered light interpolation\n \
    Measured levels of the scattered light:\n \
                           maximum {np.max(Z):.2f}\n \
                           minimum {np.min(Z):.2f}\n \
                           median {np.median(Z):.2f}\n \
                           sigma {np.std(Z):.2f}\n \
    Model of the scattered light:\n \
                           maximum {np.max(sc_light):.2f}\n \
                           minimum {np.min(sc_light):.2f}\n \
                           median {np.median(sc_light):.2f}\n \
                           sigma {np.std(sc_light):.2f}\n \
    Residual signal:\n \
                           maximum {np.max(Z_new):.2f}\n \
                           minimum  {np.min(Z_new):.2f}\n \
                           median {np.median(Z_new):.2f}\n \
                           sigma {np.std(Z_new):.2f}")


    queue.put((logging.INFO, f"File: {file_name}\n{Y.shape[1]} x {Y.shape[0]} points for scattered light interpolation\n \
    Measured levels of the scattered light:\n \
                           maximum {np.max(Z):.2f}\n \
                           minimum {np.min(Z):.2f}\n \
                           median {np.median(Z):.2f}\n \
                           sigma {np.std(Z):.2f}\n \
    Model of the scattered light:\n \
                           maximum {np.max(sc_light):.2f}\n \
                           minimum {np.min(sc_light):.2f}\n \
                           median {np.median(sc_light):.2f}\n \
                           sigma {np.std(sc_light):.2f}\n \
    Residual signal:\n \
                           maximum {np.max(Z_new):.2f}\n \
                           minimum  {np.min(Z_new):.2f}\n \
                           median {np.median(Z_new):.2f}\n \
                           sigma {np.std(Z_new):.2f}"))

    #plot
    if view:
        sc_light_med=np.median(sc_light)
        sc_light_stdv=np.std(sc_light)

        fig = matplotlib.pyplot.figure(1)
        ax = matplotlib.pyplot.gca()
        matplotlib.pyplot.imshow(sc_light, cmap=cm.Greys_r, aspect='equal',
                                     norm= matplotlib.colors.Normalize(vmin=sc_light_med-sc_light_stdv,
                                       vmax=sc_light_med+sc_light_stdv*3), interpolation='nearest')

        xxx = np.arange(0, arr_new.shape[0],1)
        yyy_new = np.median(arr_new[:, int(arr_new.shape[1]/2-10):int(arr_new.shape[1]/2+10)],1)
        yyy = np.median(arr[:, int(arr.shape[1]/2-10):int(arr.shape[1]/2+10)],1)
        fig_center_slice = matplotlib.pyplot.figure(2)
        ax_center_slice = fig_center_slice.add_subplot(111)
        matplotlib.pyplot.cla()
        ax_center_slice.plot(xxx, yyy_new, 'r')
        ax_center_slice.plot(xxx, yyy, 'b')
        matplotlib.pyplot.grid()
        matplotlib.pyplot.show()

    sc_light=np.float32(sc_light)
    hdu = fits.PrimaryHDU(sc_light)
    scatter_file = os.path.join(Path2Temp, os.path.splitext(file_name)[0] + '_SLMap.fits')
    hdu.header['SOURCE'] = file_name
    hdu.header['IMAGETYP'] = 'scattered light'
    hdu.writeto(scatter_file, overwrite=True)
    print(f"A map of the scattered light saved to {scatter_file}")
    queue.put((logging.INFO, f"A map of the scattered light saved to {scatter_file}"))

    arr_new=np.float32(arr_new)
    cleared_file = ''
    hdulist[0].data = arr_new
    cleared_file = os.path.join(dir_name, os.path.splitext(file_name)[0] + '_SLR.fits')
    if subtract:
        prihdr['HISTORY'] = 'scattered light removed'
    hdulist.writeto(cleared_file, overwrite=True)
    print(f"Cleared image has been saved to {cleared_file}")
    queue.put((logging.INFO, f"Cleared image has been saved to {cleared_file}"))
    return(cleared_file)

def process_multi(Path2Data, Path2Temp, name, ap_file, step, x_order, y_order, subtract, view, queue):
    sl_remover_data = sl_remover(Path2Data, Path2Temp, name, ap_file, step, x_order, y_order, subtract, view, queue)
    shutil.move(os.fspath(os.path.join(Path2Data, name)), os.fspath(Path2Temp))
    return sl_remover_data

def sl_subtract(Path2Data, Path2Temp, conf, subtract, queue):
    view = eval(conf['view'])
    s_flat_name = conf['s_flat_name']
    step = int(conf['step'])
    x_order = int(conf['x_order'])
    y_order = int(conf['y_order'])
    ap_file = os.path.join(Path2Temp, 'traces.txt') ##
    sl_remover_data = sl_remover(Path2Data, Path2Temp, s_flat_name, ap_file, step, x_order, y_order, subtract, view, queue)
    shutil.move(os.fspath(os.path.join(Path2Data, s_flat_name)), os.fspath(Path2Temp))
    ap_file = os.path.join(Path2Temp, 'traces.txt')

    with open(os.path.join(Path2Temp, 'obj_CRR_cleaned_list.txt'), 'r') as f:
        proc_args = [(Path2Data, Path2Temp, line.strip().split(os.sep)[-1], ap_file, step, x_order, y_order, subtract, view, queue) for line in f]
        nCPUs = os.cpu_count()
        out_list = []
        if 'threading' in conf and eval(conf['threading']) and nCPUs > 2:
            with mp.Pool(processes=nCPUs) as pool:
                res_async = pool.starmap_async(process_multi, proc_args, chunksize=nCPUs)
                res_async.wait()
                out_list.extend(res_async.get())
        else:
            for item in proc_args:
                res_mono = process_multi(Path2Data, Path2Temp, item[2], ap_file, step, x_order, y_order, subtract, view, queue)
                out_list.extend([res_mono])

    with open(os.path.join(Path2Temp, 'obj_sl_cleaned_list.txt'), 'a') as f_out:
        for item in out_list:
            print(item, file=f_out)

    flat_name = os.path.join(Path2Data, 's_flat_SLR.fits')
    list_name = os.path.join(Path2Temp, 'obj_sl_cleaned_list.txt')
    return flat_name, list_name
