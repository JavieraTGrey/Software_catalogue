from astropy.io import fits
import numpy as np

# Cut images by input and combine filters with weight


def cutting_images(line, column, filters):

    images = []
    weights = []

    for i in filters:
        im, im_wht = cut(line, column, i)
        images.append(im)
        weights.append(im_wht)

    weighted, weight_err = combined_images(images, weights)

    avgout = 'weighted.fits'
    errout = 'weighted_err.fits'

    fits.writeto(avgout, weighted, overwrite=True)
    fits.writeto(errout, weight_err, overwrite=True)

    with fits.open(filters[0] + '.fits') as source_hdulist:
        source_header = source_hdulist[0].header

    with fits.open(avgout, mode='update') as dest_hdulist:
        dest_hdulist[0].header.update(source_header)

    with fits.open(errout, mode='update') as dest_hdulist:
        dest_hdulist[0].header.update(source_header)
    return print('done')


def cut(x, y, name):
    """Receives coordinates to select a square in a FITS file.
    Input:  X: tuple of x-coordinate values
            Y: tuple of y-coordinate values
            name: file name
    Output: cropped image"""
    filename = str(name) + '.fits'
    weight_filename = str(name) + '_wht' + '.fits'
    hdu = fits.open(filename)[0]
    image = hdu.data[x[0]:x[1], y[0]:y[1]]

    weight = fits.open(weight_filename)[0]
    wht = weight.data[x[0]:x[1], y[0]:y[1]]

    return image, wht


# Combinar imagenes según el peso de cada pixel


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
    optavg = np.where(dem == 0. | np.isnan(dem), 0., num / dem)

    # Optimal error
    opterr = np.sqrt(np.where(dem <= 0 | np.isnan(dem), 0., 1. / dem))

    # Optimal noise-equalized by SNR
    comb = optavg / opterr

    return comb, opterr
