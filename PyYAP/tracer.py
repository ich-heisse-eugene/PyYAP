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
from skimage.exposure import equalize_hist

import warnings
warnings.simplefilter("ignore")
##################################################################

def order_tracer(dir_name, file_name, X_half_width, step, min_height, aperture, adaptive, view):
    poly_order = 2 # order of polynomial fit

    #read image
    hdulist = pyfits.open(dir_name+'/'+file_name)
    ordim = hdulist[0].data.copy()
    hdulist.close()

    trace_im = np.zeros((ordim.shape[0], ordim.shape[1]), dtype=np.float)

    # Search for local maxima in the cross sections and construct the mask of orders. 1 - max, 0 - rest of pixels
    for x_t in range(X_half_width, ordim.shape[1]-X_half_width, X_half_width):
        slic = np.median(ordim[:,x_t-X_half_width:x_t+X_half_width-1], 1)
        peaks_coord,_ = find_peaks(slic, height=min_height, width=2)
        trace_im[peaks_coord, x_t] = 1

    # Now rearrange the data for the futher clustering
    ones = np.array(np.where(trace_im==1))
    ord_mask = ones.T
    model = AgglomerativeClustering(n_clusters=None,
                                    affinity='euclidean',
                                    linkage='single',
                                    compute_full_tree=True,
                                    distance_threshold=5)
    model.fit(ord_mask)
    labels = model.labels_
    n_clusters = len(list(set(labels)))

    order_tab = [] # output table of traces
    width_tab = [] # output table of sizes
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
    print(f"Found {n_orders} orders")
    logging.info(f"Found {n_orders} orders")

    # Recentering orders and the final tracing
    #get points for order tracing, start from center column of image
    n_points = np.floor((ordim.shape[1]/2-X_half_width)/step)
    print(f"{1+2*n_points} points in each order for fitting")
    logging.info(f"{1+2*n_points} points in each order for fitting")
    trace_x = np.arange(ordim.shape[1]//2-X_half_width - n_points*step, ordim.shape[1]-step, step, dtype=int)

    x_coord = np.arange(ordim.shape[1])
    orders = []

    for i in range(n_orders):
        print(f"Re-trace order {i}")
        logging.info(f"Re-trace order {i}")
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
                    fwhm = 2*popt[1]*np.sqrt(2**(1/(popt[2]))-1)
                    if np.isfinite(fwhm) and fwhm > 1:
                        xfit.append(x)
                        centr.append(popt[4]+yc)
                        width.append(fwhm)
                    med_fwhm = np.median(width)
                    # print(f"Order {i}, median FWHM: {med_fwhm:.2f}, mean FWHM: {np.mean(width):.2f}")
        print(f"Fit {len(xfit)} points, median FWHM {med_fwhm:.3f}")
        logging.info(f"Fit {len(xfit)} points, median FWHM {med_fwhm:.3f}")
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

    # Output data
    width_tab = np.asarray(width_tab)
    orders = np.asarray(orders)
    text_file = open(dir_name+"/temp/traces.txt", "w") # Save results
    for i in range(orders.shape[0]):
        text_file.write("Order" + '\t' + str(i) + '\n')
        for j in x_coord:
            text_file.write(format('%.2f' % chebval(j, orders[i, :])) + '\t' + format('%.2f' % width[j]) + '\n')
    text_file.close()
    print("Data saved to " + dir_name+"/temp/traces.txt")
    logging.info("Data saved to " + dir_name+"/temp/traces.txt")

    # Display the results
    fig = plt.figure(figsize=(15, 15/(ordim.shape[1]/ordim.shape[0])), tight_layout=True)
    ax0 = fig.add_subplot(1,1,1)
    ax0.imshow(equalize_hist(ordim), cmap='gist_gray')
    ax0.set_xlabel("CCD X")
    ax0.set_ylabel("CCD Y")
    for i in range(orders.shape[0]):
        ax0.plot(x_coord, chebval(x_coord, orders[i, :]) + aperture * width_tab[i, :], 'b-', lw=0.4)
        ax0.plot(x_coord, chebval(x_coord, orders[i, :]) - aperture * width_tab[i, :], 'r-', lw=0.4)
    plt.gca().invert_yaxis()
    fig.savefig(dir_name+"/orders_map.pdf", dpi=250)
    if view:
        plt.show()

    return None

# # # ##############################test###############################
# dir_name='/home/eugene/work/mres/20201110/slit03017'
# file_name = 's_ordim.fits'
# X_half_width = 4
# min_height = 20
# view = True
# step = 10
# aperture = 1.1
# view = True
# # # ##
# order_tracer(dir_name, file_name, X_half_width, step, min_height, aperture, view)
# # ##
