import astropy.io.fits as pyfits
import numpy as np

import logging

from numpy.polynomial.chebyshev import chebfit, chebval

from scipy.optimize import curve_fit
from scipy.signal import find_peaks

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from sklearn.cluster import AgglomerativeClustering

import warnings
warnings.simplefilter("ignore")
##################################################################

def order_tracer(dir_name, file_name, X_half_width, step, min_height, aperture, adaptive, view):
    poly_order = 2 # order of polynomial fit

    #read image
    hdulist = pyfits.open(dir_name.joinpath(file_name))
    ordim = hdulist[0].data.copy()
    hdulist.close()

    trace_im = np.zeros((ordim.shape[0], ordim.shape[1]), dtype=float)

    # Search for local maxima in cross-sections and construct the mask of orders. 1 - max, 0 - rest of pixels
    for x_t in range(X_half_width, ordim.shape[1]-X_half_width, X_half_width):
        slic = np.median(ordim[:,x_t-X_half_width:x_t+X_half_width-1], 1)
        peaks_coord,_ = find_peaks(slic, height=min_height, width=2)
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
    for i in range(n_clusters):
        if len(ord_mask[labels==i,1]) >= (ordim.shape[1]/X_half_width-2)*0.5:
            cheb_coef = chebfit(ord_mask[labels==i,1], ord_mask[labels==i,0], poly_order)
            cheb_coef = np.insert(cheb_coef, 0, chebval(center_x, cheb_coef))
            order_tab.append(cheb_coef)
    order_tab = np.asarray(order_tab) # convert list into array
    order_tab = order_tab[order_tab[:,0].argsort()]  # Sort array by ascending
    order_tab[:, 0] = np.arange(len(order_tab[:, 0]), dtype=int) # insert the number of order in the 0th column
    n_orders = int(order_tab[-1, 0]) + 1
    print(f"{n_orders} orders has been detected")
    logging.info(f"{n_orders} orders has been detected")

    # Recentering orders and the final tracing
    #get points for order tracing, start from center column of image
    n_points = np.floor((ordim.shape[1]/2-X_half_width)/step)
    print(f"Each order has {1+2*n_points:.0f} points for fitting")
    logging.info(f"Each order has {1+2*n_points:.0f} points for fitting")
    trace_x = np.arange(ordim.shape[1]//2-X_half_width - n_points*step, ordim.shape[1]-step, step, dtype=int)

    x_coord = np.arange(ordim.shape[1])
    orders = []

    for i in range(n_orders):
        print(f"Re-trace order #{i}")
        logging.info(f"Re-trace order #{i}")
        xfit = []
        centr = []
        width = []
        for x in trace_x:
            xt = np.arange(x-step, x+step, 1, dtype=int)
            yc = chebval(x, order_tab[i, 1:])
            if yc > 0 and yc+X_half_width+2 < ordim.shape[0]:
                yy = np.arange(yc-X_half_width, yc+X_half_width+2, 1, dtype=int)
                prof = np.median(ordim[yy[0]:yy[-1]+1, xt], axis=1)
                if prof.shape[0] != 0 and max(prof) > 20.: # Fit only well-exposured part of the order
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
                        if np.isfinite(fwhm) and fwhm > 1 and fwhm < 10:
                            xfit.append(x)
                            centr.append(popt[4]+yc)
                            width.append(fwhm)
                        med_fwhm = np.median(width)

        print(f"{len(xfit)} points were fitted, median FWHM is {med_fwhm:.2f} pix")
        logging.info(f"{len(xfit)} points were fitted, median FWHM is {med_fwhm:.2f} pix")
        coef_center = chebfit(xfit, centr, poly_order)
        coef_width = chebfit(xfit, width, poly_order)
        ## Check the limits of the orders
        if adaptive:
            width = chebval(x_coord, coef_width)
        else:
            width = np.repeat(med_fwhm, ordim.shape[1])
        if np.min(chebval(x_coord, coef_center) - width) < 1. or np.max(chebval(x_coord, coef_center) + width) >= ordim.shape[0]-1:
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
            text_file.write(f"{i:d}\t{xcc:.3f}\t{chebval(xcc, orders[i, :]):.3f}\t{' '.join(str(jj) for jj in orders[i, :])}\t{' '.join(str(jj) for jj in width_coef_tab[i, :])}\n")
        else:
            text_file.write(f"{i:d}\t{xcc:.3f}\t{chebval(xcc, orders[i, :]):.3f}\t{' '.join(str(jj) for jj in orders[i, :])}\t{width_tab[i, 0]:.2f}\n")
    text_file.close()
    print(f"The results have been saved to {dir_name.joinpath('temp/traces.txt')}")
    logging.info(f"The results have been saved to {dir_name.joinpath('temp/traces.txt')}")

    # Display the results
    fig = plt.figure(figsize=(15, 15/(ordim.shape[1]/ordim.shape[0])), tight_layout=True)
    ax0 = fig.add_subplot(1,1,1)
    norm = matplotlib.colors.LogNorm(1, ordim.max(), clip=True)
    ax0.imshow(ordim, norm=norm, cmap='gist_gray')
    ax0.set_xlabel("CCD X")
    ax0.set_ylabel("CCD Y")
    for i in range(orders.shape[0]):
        ax0.plot(x_coord, chebval(x_coord, orders[i, :]) + aperture * width_tab[i, :], 'b-', lw=0.4)
        ax0.plot(x_coord, chebval(x_coord, orders[i, :]) - aperture * width_tab[i, :], 'r-', lw=0.4)
        ax0.text(x_coord[15], chebval(x_coord[15], orders[i, :])-3, '\#'+str(i+1), color='y', fontsize=6)
    plt.gca().invert_yaxis()
    fig.savefig(dir_name.joinpath('orders_map.pdf'), dpi=350)
    if view:
        plt.show()

    return None
