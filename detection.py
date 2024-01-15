from astropy.io import fits
from astropy.convolution import Gaussian2DKernel
import numpy as np
from scipy.optimize import curve_fit
import sep


def gaussian(x, amplitude, mean, std_dev):
    return amplitude * np.exp(-(x - mean)**2 / (2 * std_dev**2))


def bkg_stimate(data):
    # Reviso el thresh en data
    dato, bins = np.histogram(data.flatten(),
                              bins=np.linspace(-0.1, 0.005, 300))

    params, conv = curve_fit(gaussian, bins[:-1], dato)
    amplitude, mean, std_dv = params
    return amplitude, mean, std_dv


def detection(data, fwhm):
    amplitude, mean, std_dv = bkg_stimate(data)
    kernel = np.array(Gaussian2DKernel(fwhm/2.35))

    # Detection
    data = data.byteswap().newbyteorder()
    sep.set_extract_pixstack(10000000)

    objects = sep.extract(data, thresh=1.2*std_dv, minarea=3,
                          deblend_cont=0.0001, clean=False,
                          clean_param=1.0,
                          filter_kernel=kernel,
                          filter_type='matched')
    return objects, std_dv


data = fits.getdata('weighted.fits')
fwhm = 3.5
objects, std_dv = detection(data, fwhm)
