from astropy.io import fits
import astropy.units as u
from catalog_comparison import explore_square
import numpy as np
import matplotlib.pyplot as plt
import sep

names = ('weighted', 'UNCOVER_DR2_LW_SUPER_catalog', 'Images/A2744_F356W')
ceros = np.array([[(4303, 10045), (5280, 9455)], [(4477, 5095), (76, 655)]])
fwhm = 3.5
init = -5
stop = 1
thresh = 1.2
max_sep = 0.15*u.arcsec

a, b, c, d, e, f, g, h, j, k, n, o, p, q, r, t, weighted, std_dv = explore_square(names, ceros, init, stop, fwhm, thresh, max_sep)


# ordenado como matched = UNCOVER(0, 1), TOR(0, 1)
matched = c, d
unmatched_TOR = g, h
unmatched_UN = n, o


# Abro weighted
ceros = np.array([[(4303, 10045), (5280, 9455)], [(4477, 5095), (76, 655)]])

weighted_err = fits.getdata('weighted_err.fits')

photwht = weighted_err[76:655, 4477:5095]

# encuentro photerr
photerr = np.where((photwht == 0) | np.isnan(photwht),
                   np.inf, 1./np.sqrt(photwht))
photerr[~np.isfinite(photerr)] = np.median(photerr[np.isfinite(photerr)])

# recorto a cuadrito en weightes

# sum flux on circular apertures
sep.set_extract_pixstack(10000000)
weighted = weighted.copy(order='C')
weighted = weighted.byteswap().newbyteorder()


def circular_aperture(data, column, line, rad, err=photerr):
    flux, flux_err, flags = sep.sum_circle(data, column, line, rad, err)
    flux_err = np.where((flux_err == 0) | np.isnan(flux_err), 1, flux_err)
    return flux, flux_err, flags


def SNR(flux, flux_err):
    return flux / flux_err


flux_unmatch_TOR, fluxerr_unmatch_TOR, _ = circular_aperture(weighted, g, h, 2)
flux_matched, fluxerr_matched, _ = circular_aperture(weighted, c, d, 2)

SNR_TOR = SNR(flux_unmatch_TOR, fluxerr_unmatch_TOR)
SNR_matched = SNR(flux_matched, fluxerr_matched)

plt.hist(SNR_matched, bins=np.linspace(0, 100, 30), alpha=0.5)
plt.hist(SNR_TOR, bins=np.linspace(0, 100, 60), alpha=0.5)
