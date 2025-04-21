from pathlib import Path

from astropy.io import fits
import numpy as np

import logging

from scipy.optimize import curve_fit
from scipy.signal import find_peaks

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from sklearn.cluster import AgglomerativeClustering
mpl.rcParams['text.usetex'] = False

import warnings
warnings.simplefilter("ignore")
##################################################################

def order_tracer(dir_name, file_name, X_half_width, step, min_height, aperture, adaptive, view):
    poly_order = 2 # order of polynomial fit

    #read image
    hdulist = fits.open(dir_name.joinpath(file_name))
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
                                    distance_threshold=15)
    model.fit(ord_mask)
    labels = model.labels_
    n_clusters = len(list(set(labels)))

    order_tab = [] # output table of traces
    width_tab = [] # output table of sizes
    width_coef_tab = [] # output table of sizes
    center_x = ordim.shape[1]/2 - 1
    for i in range(n_clusters):
        if len(ord_mask[labels==i,1]) >= (ordim.shape[1]/X_half_width-2)*0.5:
            poly_coef = np.polyfit(ord_mask[labels==i,1], ord_mask[labels==i,0], poly_order)
            poly_coef = np.insert(poly_coef, 0, np.polyval(poly_coef, center_x))
            order_tab.append(poly_coef)
    order_tab = np.asarray(order_tab) # convert list into array
    order_tab = order_tab[order_tab[:,0].argsort()]  # Sort array by ascending
    order_tab[:, 0] = np.arange(len(order_tab[:, 0]), dtype=int) # insert the number of order in the 0th column
    n_orders = int(order_tab[-1, 0]) + 1
    print(f"{n_orders} orders has been detected")
    logging.info(f"{n_orders} orders has been detected")

    # Recentering orders and the final tracing
    #get points for order tracing, start from center column of image
    # n_points = np.floor((ordim.shape[1]/2-X_half_width)/step)
    print(f"Each order has {step} points for fitting")
    logging.info(f"Each order has {step} points for fitting")
    trace_x = np.linspace(step, ordim.shape[1]-step-1, step, dtype=int)
    x_coord = np.arange(ordim.shape[1])
    orders = []

    for i in range(n_orders):
        print(f"Re-trace order #{i}")
        logging.info(f"Re-trace order #{i}")
        xfit = []
        centr = []
        width = []
        for x in trace_x:
            xt = np.arange(x-1, x+2, 1, dtype=int)
            yc = np.polyval(order_tab[i, 1:], x)
            margin = 2.5*X_half_width
            if yc >= margin and yc < ordim.shape[0]-margin:
                yy = np.arange(yc-margin, yc+margin+1, 1, dtype=int)
                prof = np.mean(ordim[yy[0]:yy[-1]+1, xt], axis=1)
                if prof.shape[0] != 0 and max(prof) > min_height: # Fit only well-exposured part of the order
                    # Re-fit the cross-sections
                    moffat = lambda x, A, B, C, D, x0: A*(1 + ((x-x0)/B)**2)**(-C)+D
                    Y = yy-yc
                    p0 = np.array([max(prof), 3.0, 3.0, 0.0, 3.0])
                    try:
                        popt, pcov = curve_fit(moffat, Y, prof, p0, maxfev=10000)
                    except RuntimeError:
                        pass
                    else:
                        fwhm = 2*popt[1]*np.sqrt(2**(1/(popt[2]))-1)
                        if np.isfinite(fwhm) and fwhm > 1 and fwhm < 50:
                            xfit.append(x)
                            centr.append(popt[4]+yc)
                            width.append(fwhm)
                            # print(f"order = {i+1}, x = {x}, yc = {yc:.3f}, dx = {popt[4]:.3f}, fwhm = {fwhm:.3f}")
                        med_fwhm = np.median(width)

        print(f"{len(xfit)} points were fitted, median FWHM is {med_fwhm:.2f} pix")
        logging.info(f"{len(xfit)} points were fitted, median FWHM is {med_fwhm:.2f} pix")
        coef_center = np.polyfit(xfit, centr, poly_order)
        coef_width = np.polyfit(xfit, width, poly_order)
        ## Check the limits of the orders
        if adaptive:
            width = np.polyval(coef_width, x_coord)
        else:
            width = np.repeat(med_fwhm, ordim.shape[1])
        if np.min(np.polyval(coef_center, x_coord) - width) < 0. or np.max(np.polyval(coef_center, x_coord) + width) >= ordim.shape[0]-1:
            print(f"Skip incomplete order #{i}")
            logging.info(f"Skip incomplete order #{i}")
        else:
            width_tab.append(width)
            orders.append(coef_center)
            if adaptive:
                width_coef_tab.append(coef_width)

    # Output data
    width_tab = np.asarray(width_tab)
    orders = np.asarray(orders)
    text_file = open(dir_name.joinpath('temp', 'traces.txt'), "w") # Save results
    text_file.write("# N_orders   Poly_order   Adaptive?  Aperture_size_in_FWHM\n")
    text_file.write("# Order_numb   Xcent  Ycent  C(0) C(1) ... C(p) Cfwhm(0) Cfwhm(1) ... Cfwhm(p)\n")
    text_file.write("# If FWHM is fixed, Cfwhm is replaced by a single value of the median FWHM\n")
    text_file.write(f"{orders.shape[0]:d}\t{poly_order:d}\t{adaptive}\t{aperture:.1f}\n")
    for i in range(orders.shape[0]):
        xcc = np.mean(x_coord)
        if adaptive:
            text_file.write(f"{i:d}\t{xcc:.3f}\t{np.polyval(orders[i, :], xcc):.3f}\t{' '.join(str(jj) for jj in orders[i, :])}\t{' '.join(str(jj) for jj in width_coef_tab[i])}\n")
        else:
            text_file.write(f"{i:d}\t{xcc:.3f}\t{np.polyval(orders[i, :], xcc):.3f}\t{' '.join(str(jj) for jj in orders[i, :])}\t{width_tab[i, 0]:.2f}\n")
    text_file.close()
    print(f"The results have been saved to {dir_name.joinpath('temp/traces.txt')}")
    logging.info(f"The results have been saved to {dir_name.joinpath('temp/traces.txt')}")

    # Display the results
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
        ax0.text(x_coord[15], np.polyval(orders[i, :], x_coord[15])-3, '#'+str(i+1), color='y', fontsize=6)
    plt.gca().invert_yaxis()
    fig.savefig(dir_name.joinpath('orders_map.pdf'), dpi=350)
    if view:
        plt.show()

    return None
