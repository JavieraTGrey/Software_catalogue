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

line = np.array([4303, 10045])
column = np.array([5280, 9455])

# Coordenadas de las fuentes del catalogo en pixeles de la imagen original
x, y = data['x'], data['y']

# Buscar el cuadrado
select = (y < 10045) & (y > 4303) & (x > 5280) & (x < 9455)

column_real = x[select]
line_real = y[select]


# Agrupo en columna y linea
real = (line_real, column_real)

# Correccion por el pixel inicial del cuadrado
detec = (objects['y'] + line[0], objects['x'] + column[0])

# Cambio a coordenadas WCS a partir de un filtro
header = fits.getheader('Images/A2744_F356W.fits', 0)
wcs_world = WCS(header)

detcoords = wcs_world.pixel_to_world(detec[0], detec[1])
realcoords = wcs_world.pixel_to_world(real[0], real[1])


# Uso match_coordinates_sky para comparar el catalogo
index, sep2d, dist3d = match_coordinates_sky(detcoords, realcoords)
# Ordeno por Matches
catalog_matches = realcoords[index]

# ---------------------------------------------
# Visualización de las detecciones

# Primero quiero cortar la imagen weighted

cut_line = np.array([1045, 1688])
cut_column = np.array([3410, 4138])

weighted = weighted[1045:1688, 3410:4138]

# Selecciono dentro de la deteccion lo que quiero!

# Paso matches a pix
line1, column1 = wcs_world.world_to_pixel(catalog_matches)

# Cambio coordenadas de imagen recortada a imagen original

# busco el cuadrado de weighted
select1 = (line1 < line[1]) & (line1 > line[0]) & (column1 > column[0]) & (column1 < column[1])

# Cambio coordenadas al (0,0) de weighted

line1 = line1[select1] - line[0]
column1 = column1[select1] - column[0]

# Selecciono 2do recorte
select3 = (line1 < cut_line[1]) & (line1 > cut_line[0]) & (column1 < cut_column[1]) & (column1 > cut_column[0])

# Cambio coordenadas al (0,0) de 2do recorte
line2 = line1[select3] - cut_line[0]
column2 = column1[select3] - cut_column[0]


# Reviso ahora recorte en weighted

# Paso matches a pix
line3, column3 = wcs_world.world_to_pixel(detcoords)

# Cambio coordenadas de imagen recortada a imagen original

# busco el cuadrado de weighted
select4 = (line3 < line[1]) & (line3 > line[0]) & (column3 > column[0]) & (column3 < column[1])

# Cambio coordenadas al (0,0) de weighted

line3 = line3[select4] - line[0]
column3 = column3[select4] - column[0]

# Selecciono 2do recorte
select5 = (line3 < cut_line[1]) & (line3 > cut_line[0]) & (column3 < cut_column[1]) & (column3 > cut_column[0])

# Cambio coordenadas al (0,0) de 2do recorte
line4 = line3[select5] - cut_line[0]
column4 = column3[select5] - cut_column[0]

# --------

line_im = objects['y']
column_im = objects['x']

select2 = (line_im < cut_line[1]) & (line_im > cut_line[0]) & (column_im > cut_column[0]) & (column_im < cut_column[1])

# Vuelvo al origen de la imagen
line_im = line_im[select2] - cut_line[0]
column_im = column_im[select2] - cut_column[0]

# ploteo
fig, ax = plt.subplots()
ax.set_title('Todas en ambos catalogos')

m, s = np.mean(weighted), np.std(weighted)
im = ax.imshow(weighted, interpolation='nearest', cmap='gray',
               vmin=m-s, vmax=m+s, origin='lower')

ax.scatter(column2, line2, color='red', marker='*', label='UNCOVER')
ax.scatter(column4, line4, color='blue', marker='x', label='Toro')
# ax.scatter(column_im, line_im, color='yellow', marker='+', label='Check')


# ax.plot([column2, line2], [column4, line4], color='blue')
ax.legend()
# for i in range(len(objects['x']))
#     e = Ellipse(xy=(objects['x'][i], objects['y'][i]),
#                 width=6*objects['a'][i],
#                 height=6*objects['b'][i],
#                 angle=objects['theta'][i] * 180. / np.pi)
#     e.set_facecolor('none')
#     e.set_edgecolor('red')
# # plot an ellipse for each object

#     e = Ellipse(xy=(objects['x'][i], objects['y'][i]),
#                 width=6*objects['a'][i],
#                 height=6*objects['b'][i],
#                 angle=objects['theta'][i] * 180. / np.pi)
#     e.set_facecolor('none')
#     e.set_edgecolor('red')
#     ax.add_artist(e)
