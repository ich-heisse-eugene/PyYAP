from astropy.io import fits
import numpy as np
import os

from scipy.optimize import curve_fit

import logging

from scipy.optimize import curve_fit
from scipy.signal import find_peaks, correlate

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import warnings
warnings.simplefilter("ignore")

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
        ax0.plot(x_coord, np.polyval(orders[i, :], x_coord), 'g.', ms=0.5)
        ax0.plot(x_coord, np.polyval(orders[i, :], x_coord) - aperture * width_tab[i, :], 'r-', lw=0.4)
        ax0.text(x_coord[15], np.polyval(orders[i, :], x_coord[15])-3, 'Ord '+str(i+1), color='y', fontsize=6)
    plt.gca().invert_yaxis()
    fig.savefig(os.path.join(dir_name, 'orders_map.pdf'), dpi=350)
    if view:
        plt.show()
    return None

def fit_peak(img, x0, sx, y0, sy, dxy):
    x0 = x0 + sx; y0 = y0 + sy
    ind_x = np.linspace(x0-dxy, x0+dxy, 2*dxy+1, dtype=int)
    ind_y = np.linspace(y0-dxy, y0+dxy, 2*dxy+1, dtype=int)
    X, Y = np.meshgrid(ind_x, ind_y)
    chunk = img[ind_y.min():ind_y.max()+1, ind_x.min():ind_x.max()+1]
    initial_guess = (5000, x0, y0, 3, 3, 0, 0)
    try:
        popt, pcov = curve_fit(Gaussian2D, (X, Y), chunk.ravel(), p0=initial_guess)
    except Exception as e:
        print(f"Warning: {e}")
        popt = np.repeat(-99, 7)
    else:
        pass
    return popt

def Gaussian2D(xy, amp, x0, y0, sigma_x, sigma_y, theta, offset):
    """
    Adopted from https://stackoverflow.com/questions/21566379/fitting-a-2d-gaussian-function-using-scipy-optimize-curve-fit-valueerror-and-m
    """
    x, y = xy
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amp*np.exp( -(a*((x-x0)**2) + 2*b*(x-x0)*(y-y0) + c*((y-y0)**2)))
    return g.ravel()

def measure_shift(img, X0, Y0, dx, dy, dxy, view):
    plt.imshow(img, cmap="jet", origin="lower", vmin=1, vmax=1000)
    X = []; Y = []; Xnew = []; Ynew = []
    for i in range(len(X0)):
        res = fit_peak(img, X0[i], dx, Y0[i], dy, dxy)
        plt.plot(X0[i]+dx, Y0[i]+dy, 'rx', ms=3)
        if res[1] > 0 and res[2] > 0:
            X.append(X0[i]); Y.append(Y0[i])
            Xnew.append(res[1]); Ynew.append(res[2])
            plt.plot(res[1], res[2], 'y+', ms=3)
            # print(f"{res[1]:.3f}\t{res[2]:.3f}")
    plt.savefig("reflines_map.pdf", dpi=250)
    if view:
        plt.show()
    return np.asarray(X), np.asarray(Y), np.asarray(Xnew), np.asarray(Ynew)

def eval_tmatrix(x, y, xx, yy):
    S_init = np.matrix([x, y, np.ones(len(x))]).T
    S_curr = np.matrix([xx, yy, np.ones(len(x))]).T
    try:
        M = np.linalg.pinv(S_init) @ S_curr
    except Exception as e:
        print("Evaluation of the transition matrix failed")
    else:
        return M

def read_reflines(filename):
    x, y = np.loadtxt(filename, unpack=True, usecols=(0,1), comments='#')
    return x, y

def read_traces(ap_file):
    Y = []
    FWHM = []
    with open(ap_file, 'r') as fp:
        f = fp.readlines()
        p = f[3].strip().rsplit()
        if len(p) == 4:
            n_orders = int(p[0]); poly_order = int(p[1])
            adaptive = p[2]; ap_size = float(p[3])
        else:
            print("Problem of reading the file with traces")
        for i in range(4, 4+n_orders):
            p = f[i].strip().rsplit()
            poly_trace_coef = np.asarray(p[3:3+poly_order+1], dtype=float)
            poly_width_coef = np.asarray(p[3+poly_order+1:], dtype=float)
            Y.append(poly_trace_coef)
            FWHM.append(poly_width_coef)
        print(f"Read {n_orders} orders from file with traces")
    return n_orders, poly_order, adaptive, ap_size, np.asarray(Y), np.asarray(FWHM)

def rough_shift_eval(file_arc_spec, file_arc_ref):
    subar = 64
    try:
        with fits.open(file_arc_spec) as hdu:
            img_arc = hdu[0].data
    except Exception as e:
        print("Error: {e}")
        return None, None
    try:
        with fits.open(file_arc_ref) as hdu:
            img_ref = hdu[0].data
    except Exception as e:
        print("Error: {e}")
        return None, None
    (ydim, xdim) = img_ref.shape
    img_ref = img_ref[ydim//2-subar:ydim//2+subar+1, xdim//2-subar:xdim//2+subar+1]
    img_ref = img_ref / np.max(img_ref)
    img_arc = img_arc[ydim//2-subar:ydim//2+subar+1, xdim//2-subar:xdim//2+subar+1]
    img_arc = img_arc / np.max(img_arc)
    # Calculate CCF of two images
    ccf = correlate(np.float32(img_ref), np.float32(img_arc), mode="full")
    (ymax, xmax) = np.unravel_index(ccf.argmax(), ccf.shape)
    x0 = xmax - 2*subar - 1
    y0 = ymax - 2*subar - 1
    return -x0, -y0

def remap_orders(dir_name, file_arc_spec, file_reflines, file_ap, dx, dy, dxy, view, queue): # x_ref, y_ref,
    err = 1
    try:
        print(file_arc_spec)
        with fits.open(file_arc_spec) as hdu:
            img = hdu[0].data
    except Exception as e:
        print(f"Error while reading the reference image with arc lines: {e}")
        err = 0
    try:
        x_arc_ref, y_arc_ref = read_reflines(file_reflines)
    except Exception as e:
        print(f"Error while reading the coordinates of the reference lines: {e}")
        err = 0
    try:
        n_orders, poly_ord, adaptive, aperture, y_old_coef, fwhm_old_coef = read_traces(file_ap)
        y_new = np.zeros(img.shape[1])
    except Exception as e:
        print(f"Error while reading the files with apertures: {e}")
        err = 0
    coef_centr = []
    if err != 0:
        x_coord = np.arange(img.shape[1])
        X, Y, XX, YY = measure_shift(img, x_arc_ref, y_arc_ref, dx, dy, dxy, view)
        M = eval_tmatrix(X, Y, XX, YY)
        for i in range(n_orders):
            y_old = np.polyval(y_old_coef[i], x_coord)
            Vt = np.matrix([x_coord, y_old, np.ones(len(x_coord))]).T
            Vn = Vt @ M
            coef_centr.append(np.polyfit(np.asarray(Vn[:,0].T)[0], np.asarray(Vn[:,1].T)[0], poly_ord))
        coef_centr = np.asarray(coef_centr)
        ## Output data
        status = save_traces(dir_name, x_coord, poly_ord, adaptive, aperture, coef_centr, fwhm_old_coef)
        ##
        if status:
            print(f"The results have been saved to {os.path.join(dir_name, 'temp', 'traces.txt')}")
            queue.put((logging.INFO, f"The results have been saved to {os.path.join(dir_name, 'temp', 'traces.txt')}"))
        else:
            print(f"The results could not be saved to {os.path.join(dir_name, 'temp', 'traces.txt')}")
            queue.put((logging.INFO, f"The results could not be saved to {os.path.join(dir_name, 'temp', 'traces.txt')}"))
        plot_traces(dir_name, img, x_coord, coef_centr, fwhm_old_coef, aperture, view)
    return None

def remap_existing(Path2Data, Pkg_path, thars, conf, queue):
    view = eval(conf['view'])
    print("Using the existing map of orders")
    queue.put((logging.INFO, "Using the existing map of orders"))
    dx = conf['dx'] if 'dx' in conf else None
    dy = conf['dy'] if 'dy' in conf else None
    dxy = int(conf['dxy']) if 'dxy' in conf else 10
    with open(thars) as fp:
        thar_ref = fp.readline().strip('\n')
    if dx in ('None', "", None) or dy in ('None', "", None):
        dx, dy = rough_shift_eval(thar_ref, os.path.join(Pkg_path, 'devices', conf['device'], conf['device']+'_thar_ref.fits'))
    if dx != None and dy != None:
        print(f"Shift between two setups roughly estimated as dx = {dx} pix and dy = {dy} pix")
        queue.put((logging.INFO, f"Shift between two setups roughly estimated as dx = {dx} pix and dy = {dy} pix"))
        remap_orders(Path2Data, thar_ref, os.path.join(Pkg_path, 'devices', conf['device'], conf['device']+'_reflines.txt'), os.path.join(Pkg_path, 'devices', conf['device'], conf['device']+'_traces.txt'), dx, dy, dxy, view, queue)
        return True
    else:
        print("Automatic algorithm failed to evaluate shift between setups.")
        queue.put((logging.INFO, "Automatic algorithm failed to evaluate shift between setups."))
    return False
