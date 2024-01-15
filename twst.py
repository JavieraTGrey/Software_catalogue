import astropy.units as u
import numpy as np
import pyds9
from detection import explore_square
import os

names = ('weighted', 'UNCOVER_DR2_LW_SUPER_catalog', 'Images/A2744_F356W')
ceros = np.array([[(5280, 9455), (4303, 10045)], [(3410, 4138), (1045, 1688)]])
fwhm = 3.5
max_sep = 0.15*u.arcsec

a, b, c, d, e, f, g, h, j, k, n, o, p, q, r, t, weighted, std_dv = explore_square(names, ceros, fwhm, max_sep)

over = a, b, c, d
under = e, f, g, h
unmatched = j, k, n, o
all_detec = p, q, r, t

# Open a connection to DS9
ds9_path = '/home/javivi/DS9/ds9'
os.system(ds9_path + ' &')
ds9 = pyds9.DS9()

# Display the FITS image in DS9
ds9.set_np2arr(weighted)

# Define your point detections as DS9 regions (example: circles)
# Replace these with your actual coordinates and radius
point_detections = [
    {'type': 'circle', 'coords': (a, b), 'radius': 5, 'color': 'green'},
    {'type': 'circle', 'coords': (c, d), 'radius': 5, 'color': 'red'},
    # Add more detections as needed
]

# Send the DS9 regions to DS9
for detection in point_detections:
    region_str = f"{detection['type']}({detection['coords'][0]}, {detection['coords'][1]}, {detection['radius']})"
    ds9.set(region_str, color=detection['color'])

# Wait for the user to close DS9
ds9.set("zoom to fit")
ds9.set("tile")
ds9.set("raise")
