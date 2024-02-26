from astropy.io import fits
import numpy as np


def noise_equalize(filters, outname, **kwargs):
    line = kwargs.get('line', None)
    column = kwargs.get('column', None)

    images = []
    weights = []

    for i in filters:
        im, im_wht = read_data(i, x=line, y=column)
        images.append(im)
        weights.append(im_wht)  

    comb, opterr = combined_images(images, weights)

    err = str(outname) + '_opterr.fits'
    noise = str(outname) + "_noise_equal.fits"

    fits.writeto(noise, comb, overwrite=True)
    fits.writeto(err, opterr, overwrite=True)

    if line is None:
        with fits.open(filters[0] + '.fits') as source_hdulist:
            source_header = source_hdulist[0].header

        with fits.open(noise, mode='update') as dest_hdulist:
            dest_hdulist[0].header.update(source_header)

        with fits.open(err, mode='update') as dest_hdulist:
            dest_hdulist[0].header.update(source_header)
    return print('done')


def read_data(name, **kwargs):
    """Receives coordinates to select a square in a FITS file.
    Input:  X: tuple of x-coordinate values
            Y: tuple of y-coordinate values
            name: file name
    Output: cropped image"""
    
    filename = str(name) + '.fits'
    weight_filename = str(name) + '_wht' + '.fits'
    hdu = fits.open(filename)[0]
    weight = fits.open(weight_filename)[0]
    
    x = kwargs.get('x', None)
    y = kwargs.get('y', None)
    if x is None:
        image = hdu.data
        wht = weight.data
    else:
        image = hdu.data[x[0]:x[1], y[0]:y[1]]
        wht = weight.data[x[0]:x[1], y[0]:y[1]]

    return image, wht


def combined_images(images, weights):
    """Combine multiple images giving their respective pixel weights
    INPUT : images = tuple of the astronomical images
            weights = tuple of the pixel weights for each image
    OUTPUT: combined image  """

    for i in range(len(images)):
        if i == 0:
            raw_img = (images[i]) * weights[i]
            num = raw_img
            dem = weights[i]
            del raw_img
        else:
            raw_img = (images[i]) * weights[i]
            num += raw_img
            dem += weights[i]
            del raw_img

    # Optimal average
    optavg = np.where((dem == 0.) | np.isnan(dem), 0., num / dem)
    # Optimal error
    opterr = np.sqrt(np.where((dem <= 0) | (np.isnan(dem) is True), 0.,
                              1. / dem))

    # Optimal noise-equalized by SNR
    comb = optavg / opterr

    return comb, opterr