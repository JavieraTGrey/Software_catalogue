from astropy.io import fits
import astropy.units as u
from catalog_comparison import explore_square
from weight_images import cut
import numpy as np
import matplotlib.pyplot as plt
import sep

names = ('weighted', 'UNCOVER_DR2_LW_SUPER_catalog', 'Images/A2744_F356W')
ceros = np.array([[(4303, 10045), (5280, 9455)], [(4477, 5095), (76, 655)]])
fwhm = 3.5
thresh = 1.2
stop = 0.005
max_sep = 0.15*u.arcsec

a, b, c, d, e, f, g, h, j, k, n, o, p, q, r, t, weighted, std_dv = explore_square(names, ceros, stop, fwhm, thresh, max_sep)


# ordenado como matched = UNCOVER(0, 1), TOR(0, 1)
matched = a, b, c, d
unmatched_TOR = e, f, g, h
unmatched_UN = j, k, n, o
all_detec = p, q, r, t


# Abro weighted
ceros = np.array([[(4303, 10045), (5280, 9455)], [(4477, 5095), (76, 655)]])

f356w, f356w_wht = cut(ceros[0][1], ceros[0][0], 'Images/A2744_F356W')
f444w, f444w_wht = cut(ceros[0][1], ceros[0][0], 'Images/A2744_F444W')
f277w, f277w_wht = cut(ceros[0][1], ceros[0][0], 'Images/A2744_F277W')
f356w_wht = f356w_wht[76:655, 4477:5095]
f444w_wht = f444w_wht[76:655, 4477:5095]
f277w_wht = f277w_wht[76:655, 4477:5095]

# encuentro photere
weights = [f277w_wht, f356w_wht, f444w_wht]
photwht = np.sum(weights, axis=0)
photerr = np.where((photwht == 0) | np.isnan(photwht), np.inf, 1./np.sqrt(photwht))
photerr[~np.isfinite(photerr)] = np.median(photerr[np.isfinite(photerr)])

# recorto a cuadrito en weightes

# sum flux on circular apertures
sep.set_extract_pixstack(10000000)
weighted = weighted.copy(order='C')
weighted = weighted.byteswap().newbyteorder()
flux_unmatch_TOR, fluxerr_unmatch_TOR, flag = sep.sum_circle(weighted, g, h, r=2.0, err=photerr)
fluxerr_unmatch_TOR = np.where((fluxerr_unmatch_TOR == 0) | np.isnan(fluxerr_unmatch_TOR), 1, fluxerr_unmatch_TOR)

SNR_TOR = flux_unmatch_TOR / fluxerr_unmatch_TOR

plt.hist(fluxerr_unmatch_TOR, bins=np.linspace(0.018, 0.02, 100))
plt.hist(SNR_TOR, bins=np.linspace(0, 10, 100))

flux_unmatch_UN, fluxerr_unmatch_UN, flag = sep.sum_circle(weighted, r, t, r=2.0, err=photerr)

SNR_UN = flux_unmatch_UN/fluxerr_unmatch_UN

# plt.hist(fluxerr_unmatch_UN, bins=np.linspace(0.018, 0.02, 100))
plt.hist(SNR_UN, bins=np.linspace(0, 10, 100))

z, w, flag = sep.sum_circle(weighted, c, d, r=2.0, err=photerr)
SNR_matched = z/w
plt.hist(SNR_matched, bins=np.linspace(0, 100, 500))

