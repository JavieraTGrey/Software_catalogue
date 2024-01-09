from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

# Recortar imagenes por su peso


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


column = np.array([4303, 10045])
line = np.array([5280, 9455])

f277w, f277w_wht = cut(column, line, "Images/A2744_F277W")
f356w, f356w_wht = cut(column, line, "Images/A2744_F356W")
f444w, f444w_wht = cut(column, line, "Images/A2744_F444W")

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


images, weights = (f277w, f356w, f444w), (f277w_wht, f356w_wht, f444w_wht)
weighted = combined_images(images, weights)
stdv = np.std(weighted)

plt.imshow(weighted, vmin=-stdv, vmax=5*stdv)
fits.writeto('weighted.fits', weighted, overwrite=True)
fits.writeto('w.fits',
             fits.getdata('weighted_official.fits')[4303: 10045, 5280:9455],
             overwrite=True)
