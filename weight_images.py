from astropy.io import fits
import numpy as np

# Cut images by input and combine filters with weight


def cutting_images(column, line, filters):

    images = []
    weights = []

    for i in filters:
        im, im_wht = cut(column, line, i)
        images.append(im)
        weights.append(im_wht)

    weighted = combined_images(images, weights)

    fits.writeto('weighted2.fits', weighted, overwrite=True)

    with fits.open(filters[0] + '.fits') as source_hdulist:
        source_header = source_hdulist[0].header

    with fits.open('weighted2.fits', mode='update') as dest_hdulist:
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
    num = np.zeros(images[0].shape)
    dem = np.sum(weights, axis=0)
    dem[dem == 0] = 1
    for i in range(len(images)):
        num += images[i]*weights[i]
    return num / dem
