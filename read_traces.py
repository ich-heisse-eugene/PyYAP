from astropy.io import fits
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def read_traces(x_coo, ap_file):
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
        if adaptive == 'True':
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
        print(f"{n_orders} orders read from file")
    return np.asarray(Y), np.asarray(FWHM)

hdulist = fits.open(argv[1])
ordim = hdulist[0].data.copy()
hdulist.close()

aperture = 1.1

X = np.arange(ordim.shape[1])
orders, width = read_traces(X, 'traces.txt')
fig = plt.figure(figsize=(15, 15/(ordim.shape[1]/ordim.shape[0])), tight_layout=True)
ax0 = fig.add_subplot(1,1,1)
norm = matplotlib.colors.LogNorm(5, ordim.max(), clip=True)
ax0.imshow(ordim, norm=norm, cmap='gist_gray')
ax0.set_xlabel("CCD X")
ax0.set_ylabel("CCD Y")
for i in range(orders.shape[0]):
    ax0.plot(X, orders[i] + aperture * width[i], 'b-', lw=0.4)
    ax0.plot(X, orders[i] - aperture * width[i], 'r-', lw=0.4)
    ax0.text(X[15], orders[i, 15]-3, '#'+str(i+1), color='y', fontsize=6)
plt.gca().invert_yaxis()
plt.show()
