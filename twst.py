from astropy.io import fits
from astropy.convolution import Gaussian2DKernel
from astropy.coordinates.matching import match_coordinates_sky
from astropy.wcs import WCS
import numpy as np
from scipy.optimize import curve_fit
import sep
import matplotlib.pyplot as plt


weighted = fits.getdata('weighted.fits')

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
# catalog reading
cat = fits.open('UNCOVER_DR2_LW_SUPER_catalog.fits')
data = cat[1].data

# ordenado en columna, fila
real = (data['x'], data['y'])

# Paso mi deteccion a UNCOVER, (imagen original)
line = np.array([4303, 10045])
column = np.array([5280, 9455])

detec = (objects['x'] + column[0], objects['y'] + line[0])

# Convertir a ra/dec usando UNCOVER, cualquier filtro
header = fits.getheader('Images/A2744_F356W.fits', 0)
wcs_world = WCS(header)
detcoords = wcs_world.pixel_to_world(detec[0], detec[1])
realcoords = wcs_world.pixel_to_world(real[0], real[1])

# Compara deteccion entre ambos
index, sep2d, dist3d = match_coordinates_sky(detcoords, realcoords)
# reordeno realcoords
catalog_match = realcoords[index]

# Buscamos ahora el recorte
# Pasamos de ra/dec a pixeles con wcs de UNCOVER
UN_column, UN_line = wcs_world.world_to_pixel(catalog_match)
TOR_column, TOR_line = wcs_world.world_to_pixel(detcoords)

# Corregimos al cero de la imagen weighted
UN_column_w, UN_line_w = UN_column - column[0], UN_line - line[0]
TOR_column_w, TOR_line_w = TOR_column - column[0], TOR_line - line[0]

# Ahora buscamos el cuadrado seleccionado

sq_line = np.array([1045, 1688])
sq_column = np.array([3410, 4138])

select_UN = (UN_column_w > sq_column[0]) & (UN_column_w < sq_column[1]) & (UN_line_w > sq_line[0]) & (UN_line_w < sq_line[1])
select_TOR = (TOR_column_w > sq_column[0]) & (TOR_column_w < sq_column[1]) & (TOR_line_w > sq_line[0]) & (TOR_line_w < sq_line[1])

UN_column_sq, UN_line_sq = UN_column_w[select_UN], UN_line_w[select_UN]
TOR_column_sq, TOR_line_sq = TOR_column_w[select_TOR], TOR_line_w[select_TOR]

# Cambiamos a cero de cuadrado
UN_column_sq, UN_line_sq = UN_column_sq - sq_column[0], UN_line_sq - sq_line[0]
TOR_column_sq, TOR_line_sq = TOR_column_sq - sq_column[0], TOR_line_sq - sq_line[0]

# Ploteo
weighted = weighted[1045:1688, 3410:4138]
fig, ax = plt.subplots()
ax.set_title('Todas en ambos catalogos')

m, s = np.mean(weighted), np.std(weighted)
im = ax.imshow(weighted, interpolation='nearest', cmap='gray',
               vmin=m-s, vmax=m+s, origin='lower')

ax.scatter(UN_column_sq, UN_line_sq, color='red', marker='*', label='UNCOVER')
ax.scatter(TOR_column_sq, TOR_line_sq, color='blue', marker='x', label='Toro')
ax.plot([UN_column_sq[0], UN_line_sq[0]], [TOR_column_sq[0], TOR_line_sq[0]], color='green')
