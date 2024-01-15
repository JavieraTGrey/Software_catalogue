from astropy.io import fits
from astropy.convolution import Gaussian2DKernel
import numpy as np
from scipy.optimize import curve_fit
import sep


def gaussian(x, amplitude, mean, std_dev):
    ''' 1-D Gaussian function
    INPUT:  x: data
            amplitude: Function amplitude
            mean: Expected mean value
            std_dev: Value for standard deviation
    OUTPUT: 
            A*e^(-(x-mean)**2/(2*std_dev**2))
    '''
    return amplitude * np.exp(-(x - mean)**2 / (2 * std_dev**2))


def bkg_stimate(data, bkg_func, stop):
    """ Background estimation for pixel error using a provided function
    INPUT: data: Data array
           bkg_func: Function to estimate in the background.
                     Default is gaussian.
           stop: Data stop on the histogram to fit

    OUTPUT: params : fitted parameters of the bkg_func
    """
    # Reviso el thresh en data
    dato, bins = np.histogram(data.flatten(),
                              bins=np.linspace(-0.1, stop, 300))

    params, conv = curve_fit(bkg_func, bins[:-1], dato)
    return params


def detection(path, stop, fwhm, thresh, bkg_func=gaussian,
              kernel_func=Gaussian2DKernel, mask=None,
              minarea=3, deblend_count=0.0001, deblend_nthresh=32,
              clean=False, clean_param=1.0, filter_type='matched',
              segmentation_map=False):

    """ Detect sources on an fits image using sep.
    INPUT: path: path to the fits file
           stop: float to stop the data histogram to get the error on the
                 bkg estimation
           fwhm: FWHM of the kernel to use
           bkg_func(optional): Function to use for the bkg estimation.
                               Default is gaussian
           kernel_func(optional): Function to use in kernel.
                                Default is astropy.convolution.Gaussian2DKernel
           sep_params(optional): Paremeters for sep.detection
    OUTPUT:
            Objects: Extracted object parameters (structured array).
            seg(optional): Array of integers with same shape as data.
                           Only if segmentation_map is True
            std_dv: STD fitted for the bkg_func
    """
    data = fits.getdata(path)
    params = bkg_stimate(data, bkg_func, stop)
    amplitude, mean, std_dv = params
    kernel = np.array(kernel_func(fwhm/2.35))

    # BIG IMAGE
    data = data.byteswap().newbyteorder()
    sep.set_extract_pixstack(10000000)

    if segmentation_map is False:
        objects = sep.extract(data, thresh=thresh, minarea=minarea,
                              err=std_dv, mask=mask,
                              deblend_cont=deblend_count, clean=clean,
                              clean_param=clean_param,
                              filter_kernel=kernel,
                              filter_type=filter_type,
                              deblend_nthresh=deblend_nthresh,
                              segmentation_map=False)
        return objects, std_dv
    else:
        objects, seg = sep.extract(data, thresh=thresh, minarea=minarea,
                                   err=std_dv, mask=mask,
                                   deblend_cont=deblend_count, clean=clean,
                                   clean_param=clean_param,
                                   filter_kernel=kernel,
                                   filter_type=filter_type,
                                   deblend_nthresh=deblend_nthresh,
                                   segmentation_map=True)
        return objects, seg, std_dv


fwhm = 3.5
stop = 0.005
thresh = 1.2
objects, std_dv = detection('weighted.fits', stop, fwhm, thresh)
