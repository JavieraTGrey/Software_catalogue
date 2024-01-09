from astropy.io import fits
from astropy.convolution import Gaussian2DKernel
from astropy.coordinates.matching import match_coordinates_sky
from astropy.wcs import WCS
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import sep

weighted = fits.getdata('weighted.fits')
w = fits.getdata('w.fits')

# Reviso el thresh en data
data, bins = np.histogram(weighted.flatten(),
                          bins=np.linspace(-0.1, 0.005, 300))


def gaussian(x, amplitude, mean, std_dev):
    return amplitude * np.exp(-(x - mean)**2 / (2 * std_dev**2))


params, conv = curve_fit(gaussian, bins[:-1], data)
amplitude, mean, std_dv = params


kernel = np.array(Gaussian2DKernel(3.5/2.35))

# Detección
weighted = weighted.byteswap().newbyteorder()
sep.set_extract_pixstack(10000000)

objects = sep.extract(weighted, thresh=1.2*std_dv, minarea=3,
                      deblend_cont=0.0001, clean=False,
                      clean_param=1.0,
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


# catalog reading
cat = fits.open('UNCOVER_DR2_LW_SUPER_catalog.fits')
data = cat[1].data

column = np.array([4303, 10045])
line = np.array([5280, 9455])

# Coordenadas de las fuentes del catalogo en pixeles de la imagen original
x, y = data['x'], data['y']

# Buscar el cuadrado
select = (x < 10045) & (x > 4303) & (y > 5280) & (y < 9455)

line_real = x[select]
column_real = y[select]


# Agrupo en columna y linea
real = (column_real, line_real)

# Correccion por el pixel inicial del cuadrado
detec = (objects['y'] + column[0], objects['x'] + line[0])

# Cambio a coordenadas WCS a partir de un filtro 
header = fits.getheader('/Images/A2744_F356W.fits', 0)
wcs_world = WCS(header)
detcoords = wcs_world.pixel_to_world(detec[0], detec[1])
realcoords = wcs_world.pixel_to_world(real[0], real[1])


# Uso match_coordinates_sky para comparar el catalogo
index, sep2d, dist3d = match_coordinates_sky(detcoords, realcoords)
# Busco Matches
catalog_matches = realcoords[index]
