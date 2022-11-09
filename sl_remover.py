import astropy.io.fits as pyfits
import os
import shutil
import numpy as np

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
def sl_remover(dir_name, Path2Temp, file_name, ap_file, step, x_order, y_order, subtract, view):
    print(f"scatter light removing for {file_name} started")
    logging.info(f"scatter light removing for {file_name} started")
    orders = []
    Center = []

    #read apertures file
    ii=0
    with open(ap_file, 'r') as f:
        for line in f:
            one_str = line.rsplit('\t')
            if one_str[0]=='Order':
                if ii!=0:
                    Center=np.array(Center)
                    orders.append(Center)
                    Center = []
                ii=ii+1
            if one_str[0]!='Order':
                Center.append(float(one_str[0]))

    Center=np.asarray(Center)
    orders.append(Center)
    f.close()

    print(f"{ii} orders read from file")

    orders = np.asarray(orders)

    X = []
    #calculate middle line for 2 close orders
    for jj in range(0, orders.shape[0]-1):
        X.append((orders[jj, :]+orders[jj+1, :])/2)
    X = np.asarray(X)   #array with gaps

    background_X = np.arange(step/2, X.shape[1], step) #array with Xcoo of background probe[number of probe along X]

    Y = []
    for ii in range(0, background_X.shape[0]):
        Y.append(X[:, int(background_X[ii])])

    Y = np.asarray(Y).T #array with Ycoo of background probe[number of probe along X, number of probe along Y]

    print(f"{Y.shape[1]} x {Y.shape[0]} points for scatter light interpolation")
    logging.info(f"{Y.shape[1]} x {Y.shape[0]} points for scatter light interpolation")

    #read image file
    hdulist = pyfits.open(dir_name.joinpath(file_name))
    arr = hdulist[0].data.copy()
    prihdr = hdulist[0].header
    hdulist.close()

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
    print(f"Measured scatter light:")
    print(f"maximum {np.max(Z):.2f}")
    print(f"minimum {np.min(Z):.2f}")
    print(f"median {np.median(Z):.2f}")
    print(f"sigma {np.std(Z)}")
    logging.info(f"Measured scatter light:")
    logging.info(f"maximum {np.max(Z):.2f}")
    logging.info(f"minimum {np.min(Z):.2f}")
    logging.info(f"median {np.median(Z):.2f}")
    logging.info(f"sigma {np.std(Z):.2f}")
    print()

    # # Debug. Surface plot of scattered light
    if view:
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.plot_surface(XX, Y, Z, cmap='plasma', edgecolor='none')
        ax.set_title('Surface plot')
        plt.show()

    # sc_light = median_filter(surf_fit(XX, Y, Z, x_order, y_order, arr.shape[0], arr.shape[1]), size=15)
    sc_light = surf_fit(XX, Y, Z, x_order, y_order, arr.shape[0], arr.shape[1])

    print("Model scatter light:")
    print(f"maximum {np.max(sc_light):.2f}")
    print(f"minimum {np.min(sc_light):.2f}")
    print(f"median {np.median(sc_light):.2f}")
    print(f"sigma {np.std(sc_light):.2f}")
    logging.info("Model scatter light:")
    logging.info(f"maximum {np.max(sc_light):.2f}")
    logging.info(f"minimum {np.min(sc_light):.2f}")
    logging.info(f"median {np.median(sc_light):.2f}")
    logging.info(f"sigma {np.std(sc_light):.2f}")
    print()

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

    print("Residual scatter light:")
    print(f"maximum {np.max(Z_new):.2f}")
    print(f"minimum  {np.min(Z_new):.2f}")
    print(f"median {np.median(Z_new):.2f}")
    print(f"sigma {np.std(Z_new):.2f}")
    logging.info(f"maximum {np.max(Z_new):.2f}")
    logging.info(f"minimum  {np.min(Z_new):.2f}")
    logging.info(f"median {np.median(Z_new):.2f}")
    logging.info(f"sigma {np.std(Z_new):.2f}")
    print()

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
    hdu = pyfits.PrimaryHDU(sc_light)
    scatter_file = Path2Temp.joinpath(os.path.splitext(file_name)[0] + '_SLMap.fits')
    hdu.header['SOURCE'] = file_name
    hdu.header['IMAGETYP'] = 'scattered light'
    hdu.writeto(scatter_file, overwrite=True)
    print(f"Map of the scattered light was saved to {scatter_file}")
    logging.info(f"Map of the scattered light was saved to {scatter_file}")

    arr_new=np.float32(arr_new)
    cleared_file = ''
    hdulist[0].data = arr_new
    cleared_file = dir_name.joinpath(os.path.splitext(file_name)[0] + '_SLR.fits')
    if subtract:
        prihdr['HISTORY'] = 'scatter light removed'
    hdulist.writeto(cleared_file, overwrite=True)
    print(f"cleared image saved to {cleared_file}")
    logging.info(f"cleared image saved to {cleared_file}")

    return(cleared_file)
