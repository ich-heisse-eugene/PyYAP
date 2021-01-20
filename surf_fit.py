import numpy as np
from numpy.polynomial.chebyshev import chebfit, chebval
from scipy.interpolate import splrep, splev

import matplotlib.pyplot as plt

def fit_y(X, Y, Z, y_order, y_shape):
    y_model = []
    y_coo = np.arange(0, y_shape, 1)
    for ii in range(0, X.shape[1]):
        # tck = chebfit(Y[:, ii], Z[:, ii], y_order)
        # zy_model = chebval(y_coo, tck)
        tck = splrep(Y[:, ii], Z[:, ii], k=y_order)
        zy_model = splev(y_coo, tck)
        y_model.append(zy_model)
    return np.array(y_model)

def fit_x(y_model, x4y, x_shape, x_order):
    sc_light = []
    x_coo = np.arange(0, x_shape, 1)
    for ii in range(0, y_model.shape[1]):
        # tck = chebfit(x4y, y_model[:, ii], x_order)
        # zx_model = chebval(x_coo, tck)
        tck = splrep(x4y, y_model[:, ii], k=x_order)
        zx_model = splev(x_coo, tck)
        sc_light.append(zx_model)
    return np.array(sc_light)

def surf_fit(X, Y, Z, x_order, y_order, y_shape, x_shape):
    y_model = fit_y(X, Y, Z, y_order, y_shape)
    sc_light = fit_x(y_model, X[0], x_shape, x_order)
    return sc_light
