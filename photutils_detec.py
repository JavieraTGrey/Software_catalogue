from detection import bkg_estimate
from astropy.convolution import convolve
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
import numpy as np
from photutils.segmentation import SourceFinder, SourceCatalog
import matplotlib.pyplot as plt
from astropy.convolution import Gaussian2DKernel
from catalog_comparison import offset, catalog_comparison, get_square_on_image
import pyds9

# input data
weighted_image = fits.open('weighted.fits')
weighted = weighted_image[0].data
ceros = np.array([[(4303, 10045), (5280, 9455)], [(4477, 5103), (76, 668)]])
line, column = ceros[1]
data = weighted[column[0]:column[1], line[0]:line[1]]

UNCOVER = fits.open('UNCOVER_DR2_LW_SUPER_catalog.fits')
UNCOVER_data = UNCOVER[1].data
real = (UNCOVER_data['x'], UNCOVER_data['y'])

# PARAMS
fwhm = 3.5
init = -5.
stop = 0.8
thresh = 1.2
params = bkg_estimate(data, init, stop)
std_dv = params[2]
max_sep = 0.2*u.arcsec  # aumente un poco


# Kernel
kernel = np.array(Gaussian2DKernel(fwhm/2.35))
convolved_data = convolve(data, kernel)

# Segment map check
deblend_thresh = 32
deblend_count = 0.001
minarea = 3
finder = SourceFinder(npixels=minarea, deblend=True, nlevels=32,
                      contrast=deblend_count, progress_bar=True)
segm = finder(convolved_data, thresh*std_dv)

# fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12.5))
# ax1.imshow(data, origin='lower', cmap='Greys_r',
#            vmin=-std_dv, vmax=15*std_dv)
# ax1.set_title('Background-subtracted Data')
# ax2.imshow(segm.data, origin='lower', cmap=segm.cmap,
#            interpolation='nearest')
# ax2.set_title('Segmentation Image')


# Catalog
# AQUI DATA DEBERIA SER PSF MATCHED PERO LO HARÉ DESPUES
cat = SourceCatalog(data, segm, convolved_data=convolved_data,
                    progress_bar=True)
objects = cat.xcentroid, cat.ycentroid

wcs_world = WCS(weighted_image[0].header)
detec = offset(ceros[1], objects, '+')
detec = offset(ceros[0], detec, '+')

detcoords = wcs_world.pixel_to_world(detec[0], detec[1])
realcoords = wcs_world.pixel_to_world(real[0], real[1])

catalog_match, under, over, unmatched = catalog_comparison(max_sep,
                                                           detcoords,
                                                           realcoords)
match_under, detcoords_under = under
match_over,  detcoords_over = over
unmatch_real, unmatch_det = unmatched

# catalog_0, catalog_1, matched_0, matched_1
a, b, c, d = get_square_on_image(detcoords_under, match_under, ceros[0],
                                 ceros[1], wcs_world)
e, f, g, h = get_square_on_image(detcoords_over, match_over, ceros[0],
                                 ceros[1], wcs_world)
j, k, n, o = get_square_on_image(detcoords, catalog_match, ceros[0],
                                 ceros[1], wcs_world)
p, q, r, t = get_square_on_image(unmatch_real, unmatch_det, ceros[0],
                                 ceros[1], wcs_world)

# Open a connection to DS9
ds9 = pyds9.DS9()
ds9.set_np2arr(data)
ds9.set("scale log")
ds9.set("zoom to fit")

# Perfect detection
for i in range(len(a)):
    ds9.set("region command {cross point " + str(a[i]) + " " + str(b[i]) + " }")
    ds9.set('region select all')
    ds9.set('region color Red')

for i in range(len(c)):
    ds9.set("region command {circle " + str(c[i]) + " " + str(d[i]) + " 2}")


# Toro unmatched
ds9.set('frame new')
ds9.set_np2arr(data)
ds9.set("scale log")
ds9.set("zoom to fit")
for i in range(len(g)):
    ds9.set("region command {circle " + str(g[i]) + " " + str(h[i]) + " 2}")

# UNCOVER unmatched
ds9.set('frame new')
ds9.set_np2arr(data)
ds9.set("scale log")
ds9.set("zoom to fit")
for i in range(len(p)):
    ds9.set("region command {cross point " + str(p[i]) + " " + str(q[i]) + " }")
    ds9.set('region select all')
    ds9.set('region color Red')

# ALL
ds9.set('frame new')
ds9.set_np2arr(data)
ds9.set("scale log")
ds9.set("zoom to fit")

for i in range(len(j)):
    ds9.set("region command {cross point " + str(j[i]) + " " + str(k[i]) + " }")
    ds9.set('region select all')
    ds9.set('region color Red')
for i in range(len(p)):
    ds9.set("region command {cross point " + str(p[i]) + " " + str(q[i]) + " }")
    ds9.set('region select all')
    ds9.set('region color Red')

for i in range(len(n)):
    ds9.set("region command {circle " + str(n[i]) + " " + str(o[i]) + " 2}")

ds9.set('frame lock physical')
ds9.set('region select none')
