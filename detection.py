from astropy.io import fits
from astropy.convolution import Gaussian2DKernel, Tophat2DKernel
import numpy as np
import sep
from scipy.ndimage import gaussian_filter
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt


weighted = fits.getdata('weighted.fits')
w = fits.getdata('w.fits')

# Hago una detección de las fuentes

# Creación de kernel gaussiano
fwhm = 3.5  # 3.5 pixeles
sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
gaussian_kernel = gaussian_filter(np.ones((10, 10)), sigma)

# Normalizo kernel
gaussian_kernel /= np.sum(gaussian_kernel)

# detección
weighted = weighted.byteswap().newbyteorder()
sep.set_extract_pixstack(weighted.size)

# thresh debe ser del orden del data!

objects = sep.extract(weighted, thresh=1e-9, minarea=3,
                      deblend_cont=0.0001, clean='N',
                      filter_kernel=gaussian_kernel,
                      filter_type='matched')


# Visualización de fuentes detectadas
fig, ax = plt.subplots()
m, s = np.mean(weighted), np.std(weighted)
im = ax.imshow(weighted, interpolation='nearest', cmap='gray',
               vmin=-1e-09, vmax=4e-8, origin='lower')

# plot an ellipse for each object
for i in range(len(objects)):
    e = Ellipse(xy=(objects['x'][i], objects['y'][i]),
                width=6*objects['a'][i],
                height=6*objects['b'][i],
                angle=objects['theta'][i] * 180. / np.pi)
    e.set_facecolor('none')
    e.set_edgecolor('red')
    ax.add_artist(e)
