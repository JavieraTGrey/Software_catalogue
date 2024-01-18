import astropy.units as u
import numpy as np
import pyds9
from catalog_comparison import explore_square

names = ('weighted', 'UNCOVER_DR2_LW_SUPER_catalog', 'Images/A2744_F356W')
ceros = np.array([[(4303, 10045), (5280, 9455)], [(4477, 5095), (76, 655)]])
fwhm = 3.5
init = -5
stop = 1
thresh = 1.2
max_sep = 0.15*u.arcsec

a, b, c, d, e, f, g, h, j, k, n, o, p, q, r, t, weighted, std_dv = explore_square(names, ceros, init, stop, fwhm, thresh, max_sep)

over = a, b, c, d
under = e, f, g, h
unmatched = j, k, n, o
all_detec = p, q, r, t

# Open a connection to DS9
ds9 = pyds9.DS9()
ds9.set_np2arr(weighted)
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
ds9.set_np2arr(weighted)
ds9.set("scale log")
ds9.set("zoom to fit")
for i in range(len(g)):
    ds9.set("region command {circle " + str(g[i]) + " " + str(h[i]) + " 2}")

# UNCOVER unmatched
ds9.set('frame new')
ds9.set_np2arr(weighted)
ds9.set("scale log")
ds9.set("zoom to fit")
for i in range(len(r)):
    ds9.set("region command {cross point " + str(r[i]) + " " + str(t[i]) + " }")
    ds9.set('region select all')
    ds9.set('region color Red')

# ALL
ds9.set('frame new')
ds9.set_np2arr(weighted)
ds9.set("scale log")
ds9.set("zoom to fit")

for i in range(len(j)):
    ds9.set("region command {cross point " + str(j[i]) + " " + str(k[i]) + " }")
    ds9.set('region select all')
    ds9.set('region color Red')
for i in range(len(r)):
    ds9.set("region command {cross point " + str(r[i]) + " " + str(t[i]) + " }")
    ds9.set('region select all')
    ds9.set('region color Red')

for i in range(len(n)):
    ds9.set("region command {circle " + str(n[i]) + " " + str(o[i]) + " 2}")

ds9.set('frame lock physical')
ds9.set('region select none')
