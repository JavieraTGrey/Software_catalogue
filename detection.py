from astropy.io import fits
from astropy.convolution import Gaussian2DKernel
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import norm
import sep

weighted = fits.getdata('weighted.fits')
w = fits.getdata('w.fits')

# Reviso el thresh en data
data, bins = np.histogram(weighted.flatten(),
                          bins=np.linspace(-0.1, 0.005, 300))


def gaussian(x, amplitude, mean, std_dev):
    return amplitude * np.exp(-(x - mean)**2 / (2 * std_dev**2))


params, conv = curve_fit(gaussian, bins[:-1], data)
amplitude, mean, std = params
thresh = 2 * np.sqrt(2 * np.log(2)) * std

kernel = np.array(Gaussian2DKernel(3.5/2.35))

# detección
weighted = weighted.byteswap().newbyteorder()
sep.set_extract_pixstack(10000000)

objects = sep.extract(weighted, thresh=thresh, minarea=3,
                      deblend_cont=0.0001, clean=False,
                      clean_param=1.0,
                      err=std,
                      filter_kernel=kernel,
                      filter_type='matched')


# Visualización de las detecciones
# fig, ax = plt.subplots()
# m, s = np.mean(weighted), np.std(weighted)
# im = ax.imshow(weighted, interpolation='nearest', cmap='gray',
#                vmin=m-s, vmax=m+s, origin='lower')
# # plot an ellipse for each object
# for i in range(len(objects)):
#     e = Ellipse(xy=(objects['x'][i], objects['y'][i]),
#                 width=6*objects['a'][i],
#                 height=6*objects['b'][i],
#                 angle=objects['theta'][i] * 180. / np.pi)
#     e.set_facecolor('none')
#     e.set_edgecolor('red')
#     ax.add_artist(e)
