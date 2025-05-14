import os
import shutil

from astropy.io import fits
import numpy as np

import logging

from scipy.optimize import curve_fit
from scipy.signal import find_peaks

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from medianer import medianer

from sklearn.cluster import AgglomerativeClustering

import warnings
warnings.simplefilter("ignore")
##################################################################

def save_traces(dir_name, x_coord, poly_order, adaptive, aperture, orders, width_coef_tab):
    text_file = open(os.path.join(dir_name, 'temp', 'traces.txt'), "w") # Save results
    text_file.write("# N_orders   Poly_order   Adaptive?  Aperture_size_in_FWHM\n")
    text_file.write("# Order_numb   Xcent  Ycent  C(0) C(1) ... C(p) Cfwhm(0) Cfwhm(1) ... Cfwhm(p)\n")
    text_file.write("# If FWHM is fixed, Cfwhm is replaced by a single value of the median FWHM\n")
    text_file.write(f"{orders.shape[0]:d}\t{poly_order:d}\t{adaptive}\t{aperture:.1f}\n")
    for i in range(orders.shape[0]):
        xcc = np.mean(x_coord)
        if adaptive:
            text_file.write(f"{i:d}\t{xcc:.3f}\t{np.polyval(orders[i, :], xcc):.3f}\t{' '.join(str(jj) for jj in orders[i, :])}\t{' '.join(str(jj) for jj in width_coef_tab[i])}\n")
        else:
            text_file.write(f"{i:d}\t{xcc:.3f}\t{np.polyval(orders[i, :], xcc):.3f}\t{' '.join(str(jj) for jj in orders[i, :])}\t{width_coef_tab[i]:.2f}\n")
    text_file.close()
    return os.path.isfile(os.path.join(dir_name, 'temp', 'traces.txt'))

def plot_traces(dir_name, ordim, x_coord, orders, width_tab, aperture, view):
    fig = plt.figure(figsize=(15, 15/(ordim.shape[1]/ordim.shape[0])), tight_layout=True)
    ax0 = fig.add_subplot(1,1,1)
    norm = mpl.colors.LogNorm(1, ordim.max(), clip=True)
    ax0.imshow(ordim, norm=norm, cmap='gist_gray')
    ax0.set_xlabel("CCD X")
    ax0.set_ylabel("CCD Y")
    for i in range(orders.shape[0]):
        ax0.plot(x_coord, np.polyval(orders[i, :], x_coord) + aperture * width_tab[i, :], 'b-', lw=0.4)
        ax0.plot(x_coord, np.polyval(orders[i, :], x_coord), 'gx', ms=1)
        ax0.plot(x_coord, np.polyval(orders[i, :], x_coord) - aperture * width_tab[i, :], 'r-', lw=0.4)
        ax0.text(x_coord[15], np.polyval(orders[i, :], x_coord[15])-3, 'Ord '+str(i+1), color='y', fontsize=6)
    plt.gca().invert_yaxis()
    fig.savefig(os.path.join(dir_name, 'orders_map.pdf'), dpi=350)
    if view:
        plt.show()
    return None


def order_tracer(dir_name, file_name, X_half_width, step, min_height, aperture, adaptive, view, queue):
    poly_order = 2 # order of polynomial fit

    #read image
    with fits.open(os.path.join(dir_name, file_name)) as hdulist:
        ordim = hdulist[0].data.copy()
        hdulist.close()

    trace_im = np.zeros((ordim.shape[0], ordim.shape[1]), dtype=float)

    # Search for local maxima in cross-sections and construct the mask of orders. 1 - max, 0 - rest of pixels
    for x_t in range(3, ordim.shape[1]-3, 3):
        med_cs = np.mean(ordim[:,x_t-3:x_t+4], 1)
        slic = med_cs / np.max(med_cs) * 100.
        peaks_coord,_ = find_peaks(slic, height=min_height, distance=1.22*X_half_width)
        trace_im[peaks_coord, x_t] = 1

    # Now re-arrange the data for futher clustering
    ones = np.array(np.where(trace_im==1))
    ord_mask = ones.T
    model = AgglomerativeClustering(n_clusters=None,
                                    metric='euclidean',
                                    linkage='single',
                                    compute_full_tree=True,
                                    distance_threshold=5)
    model.fit(ord_mask)
    labels = model.labels_
    n_clusters = len(list(set(labels)))

    order_tab = [] # output table of traces
    width_tab = [] # output table of sizes
    width_coef_tab = [] # output table of sizes
    center_x = ordim.shape[1]/2 - 1
    # plt.imshow(ordim, cmap="jet", origin="lower")
    for i in range(n_clusters):
        if len(ord_mask[labels==i,1]) >= (ordim.shape[1]/X_half_width-2)*0.5:
            # plt.plot(ord_mask[labels==i,1], ord_mask[labels==i,0], marker='.', ms=2, ls='')
            poly_coef = np.polyfit(ord_mask[labels==i,1], ord_mask[labels==i,0], poly_order)
            poly_coef = np.insert(poly_coef, 0, np.polyval(poly_coef, center_x))
            # plt.plot(np.arange(ordim.shape[1]), np.polyval(poly_coef[1:], np.arange(ordim.shape[1])), 'y--', lw=1)
            order_tab.append(poly_coef)
    # plt.show()
    order_tab = np.asarray(order_tab) # convert list into array
    order_tab = order_tab[order_tab[:,0].argsort()]  # Sort array by ascending
    order_tab[:, 0] = np.arange(len(order_tab[:, 0]), dtype=int) # insert the number of order in the 0th column
    n_orders = int(order_tab[-1, 0]) + 1
    print(f"{n_orders} orders have been detected")
    queue.put((logging.INFO, f"{n_orders} orders have been detected"))

    # Recentering orders and the final tracing
    # get points for order tracing, start from center column of image
    # n_points = np.floor((ordim.shape[1]/2-X_half_width)/step)
    trace_x = np.arange(step, ordim.shape[1]-step+1, step, dtype=int)
    n_trace_x = len(trace_x)
    print(f"Each order has {n_trace_x} points for fitting")
    queue.put((logging.INFO, f"Each order has {n_trace_x} points for fitting"))
    x_coord = np.arange(ordim.shape[1])
    orders = []

    for i in range(n_orders):
        print(f"Re-trace order #{i}")
        queue.put((logging.INFO, f"Re-trace order #{i}"))
        xfit = []
        centr = []
        width = []
        med_fwhm = 0
        for x in trace_x:
            xt = np.arange(x-1, x+2, 1, dtype=int)
            yc = np.polyval(order_tab[i, 1:], x)
            margin = X_half_width
            if yc >= margin and yc < ordim.shape[0]-margin:
                yy = np.arange(yc-margin, yc+margin+1, 1, dtype=int)
                prof = np.mean(ordim[yy[0]:yy[-1]+1, xt], axis=1); prof = prof - np.min(prof)
                prof = prof / np.max(prof)
                Y = yy-yc
                prof_at_zero = np.mean(prof[int(np.min(Y))-1:int(np.min(Y))+2])
                if prof.shape[0] > 4:
                    # Re-fit the cross-sections
                    moffat = lambda x, A, B, C, D, x0: A*(1 + ((x-x0)/B)**2)**(-C)+D
                    p0 = np.array([prof_at_zero, 3.0, 3.0, 0.0, 3.0])
                    try:
                        popt, pcov = curve_fit(moffat, Y, prof, p0, maxfev=10000)
                    except RuntimeError:
                        pass
                    else:
                        fwhm = 2*popt[1]*np.sqrt(2**(1/(popt[2]))-1)
                        if np.isfinite(fwhm) and fwhm > 1 and np.abs(popt[4]) <= X_half_width/2:
                            xfit.append(x)
                            centr.append(popt[4]+yc)
                            width.append(fwhm)
                            # plt.plot(Y, prof, ls='-', marker='s', ms=3)
                            # plt.plot(popt[4], 0, 'kX', ms=5)
                        med_fwhm = np.median(width)
        # plt.show()

        print(f"{len(xfit)} points were fitted, median FWHM is {med_fwhm:.2f} pix")
        queue.put((logging.INFO, f"{len(xfit)} points were fitted, median FWHM is {med_fwhm:.2f} pix"))
        if len(xfit) > 0:
            coef_center = np.polyfit(xfit, centr, poly_order)
            coef_width = np.polyfit(xfit, width, poly_order)
            ## Check the limits of the orders
            if adaptive:
                width = np.polyval(coef_width, x_coord)
            else:
                width = np.repeat(med_fwhm, ordim.shape[1])
            if np.min(np.polyval(coef_center, x_coord) - width) < 0. or np.max(np.polyval(coef_center, x_coord) + width) >= ordim.shape[0]-1 or len(xfit) < 0.66*n_trace_x:
                print(f"Skip incomplete order #{i}")
                queue.put((logging.INFO, f"Skip incomplete order #{i}"))
            else:
                width_tab.append(width)
                orders.append(coef_center)
                if adaptive:
                    width_coef_tab.append(coef_width)
                else:
                    width_coef_tab.append(med_fwhm)
        else:
            print(f"Skip incomplete order #{i}")
            queue.put((logging.INFO, f"Skip incomplete order #{i}"))

    width_tab = np.asarray(width_tab)
    orders = np.asarray(orders)
    width_coef_tab = np.asarray(width_coef_tab)
    ## Output data
    status = save_traces(dir_name, x_coord, poly_order, adaptive, aperture, orders, width_coef_tab)
    ##
    if status:
        print(f"The results have been saved to {os.path.join(dir_name, 'temp', 'traces.txt')}")
        queue.put((logging.INFO, f"The results have been saved to {os.path.join(dir_name, 'temp', 'traces.txt')}"))
    else:
        print(f"The results could not be saved to {os.path.join(dir_name, 'temp', 'traces.txt')}")
        queue.put((logging.INFO, f"The results could not be saved to {os.path.join(dir_name, 'temp', 'traces.txt')}"))

    # Display the results
    plot_traces(dir_name, ordim, x_coord, orders, width_tab, aperture, view)
    return None

def locate_orders(Path2Data, Path2Temp, conf, s_ordim_method, queue):
    slice_half_width = float(conf['slice_half_width'])
    step = int(conf['step'])
    min_height = float(conf['min_height'])
    aperture = float(conf['aperture'])
    view = eval(conf['view'])
    adaptive = eval(conf['adaptive'])
    s_flat_name = conf['s_flat_name']
    if 's_ordim_name' in conf:
        s_ordim_name = conf['s_ordim_name'].rstrip()
    else:
        s_ordim_name = 's_ordim.fits'
    objects_list = []
    if s_ordim_method == "hybrid" or s_ordim_method == "objects":
        with open(os.path.join(Path2Temp, 'obj_CRR_cleaned_list.txt'), 'r') as ff:
            objects_list = ff.read().splitlines()
            ff.close()
        if s_ordim_method == "hybrid":
            objects_list.append(os.path.join(Path2Data, s_flat_name))
        with open(os.path.join(Path2Temp, 'ordim_list.txt'), 'w+') as f:
            print(*objects_list, sep='\n', file=f)
            f.close()
        sordim_data = medianer(Path2Data, os.path.join(Path2Temp, 'ordim_list.txt'), os.path.join(Path2Data, s_ordim_name))
    if s_ordim_method == 'flats':
        shutil.copy2(os.path.join(Path2Data, s_flat_name), os.path.join(Path2Data, s_ordim_name))
    if os.path.isfile(os.path.join(Path2Data, s_ordim_name)):
        print(f"Master image {s_ordim_name} for the tracer was created using the method '{s_ordim_method}'")
        queue.put((logging.INFO, f"Master image {s_ordim_name} for the tracer was created using the method '{s_ordim_method}'"))
    else:
        print(f"Error: Failed to create a file {s_ordim_name} for the tracing using the method '{s_ordim_method}'")
        queue.put((logging.INFO, f"Error: Failed to create a file {s_ordim_name} for the tracing using the method '{s_ordim_method}'"))
        return False
    ## Trace orders from scratch using cluster analysis
    print("Start tracing orders")
    order_tracer(Path2Data, s_ordim_name, slice_half_width, step, min_height, aperture, adaptive, view, queue)
    shutil.move(os.fspath(os.path.join(Path2Data, s_ordim_name)), os.fspath(Path2Temp))
    return True
