from astropy.io import fits
from astropy.convolution import Gaussian2DKernel
import numpy as np
from scipy.optimize import curve_fit
import sep


def gaussian(x, amplitude, mean, std_dev):
    ''' 1-D Gaussian function.

    INPUT:
    - x (array): Data.
    - amplitude (float): Function amplitude.
    - mean (float): Expected mean value.
    - std_dev (float): Standard deviation.

    OUTPUT:
    - array: A * exp(-(x - mean)**2 / (2 * std_dev**2))
    '''
    return amplitude * np.exp(-(x - mean)**2 / (2 * std_dev**2))


def bkg_estimate(data, bkg_func=None, stop=0.005):
    """Estimate background using a provided function.

    INPUT:
    - data (array): Data array.
    - bkg_func (function, optional): Background estimation function.
                                     Default is Gaussian.
    - stop (float, optional): Data stop on the histogram to fit.
                              Default is 0.005.

    OUTPUT:
    - params (array): Fitted parameters of the background function.
    """
    if bkg_func is None:
        bkg_func = gaussian  # Default background function

    # Calculate histogram
    hist_data, bins = np.histogram(data.flatten(),
                                   bins=np.linspace(-0.1, stop, 300))

    # Fit background function to the histogram
    params, _ = curve_fit(bkg_func, bins[:-1], hist_data)
    return params


def detection(path, stop, fwhm, thresh,
              bkg_func=gaussian, err=None,
              kernel_func=Gaussian2DKernel, mask=None,
              minarea=3, deblend_count=0.0001, deblend_nthresh=32,
              clean=False, clean_param=1.0, filter_type='matched',
              segmentation_map=False):

    """Detect sources on a fits image using sep.

    INPUT:
        - path (str): Path to the fits file.
        - stop (float): Float to stop the data histogram
                        for estimating the background.
        - fwhm (float): FWHM of the kernel to use.
        - bkg_func (function, optional): Function for background estimation.
                                         Default is gaussian.
        - kernel_func (function, optional): Function for the kernel.
                    Default is astropy.convolution.Gaussian2DKernel.
        - err (array-like, optional): Error array. Default is None.
        - mask (array-like, optional): Mask array. Default is None.
        - minarea (int, optional): Minimum area for detected objects.
                                   Default is 3.
        - deblend_count (float, optional): Deblend constant. Default is 0.0001.
        - deblend_nthresh (int, optional): Deblend threshold. Default is 32.
        - clean (bool, optional): Perform cleaning. Default is False.
        - clean_param (float, optional): Cleaning parameter. Default is 1.0.
        - filter_type (str, optional): Filter type. Default is 'matched'.
        - segmentation_map (bool, optional): Return segmentation map.
                                             Default is False.

    OUTPUT:
        - objects (structured array): Extracted object parameters.
        - seg (array): Segmentation map (if segmentation_map is True).
        - std_dv (float): STD fitted for the bkg_func.
    """
    try:
        data = fits.getdata(path)
    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {path}")
    kernel = np.array(kernel_func(fwhm/2.35))

    # BIG IMAGE
    data = data.byteswap().newbyteorder()
    sep.set_extract_pixstack(10000000)
    if err is None:
        params = bkg_estimate(data, bkg_func, stop)
        amplitude, mean, std_dv = params
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
    else:
        if segmentation_map is False:
            objects = sep.extract(data, thresh=thresh, minarea=minarea,
                                  err=err, mask=mask,
                                  deblend_cont=deblend_count, clean=clean,
                                  clean_param=clean_param,
                                  filter_kernel=kernel,
                                  filter_type=filter_type,
                                  deblend_nthresh=deblend_nthresh,
                                  segmentation_map=False)
            return objects, std_dv
        else:
            objects, seg = sep.extract(data, thresh=thresh, minarea=minarea,
                                       err=err, mask=mask,
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

# ----------
# Trying with CEERS DATA
