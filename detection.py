# Hago una detección de las fuentes

# Creación de kernel gaussiano
fwhm = 3.5  # 3.5 pixeles
sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
gaussian_kernel = gaussian_filter(np.ones((10, 10)), sigma)

# Normalizo kernel
gaussian_kernel /= np.sum(gaussian_kernel)

# detección
objects, segmap = sep.extract(weighted, thresh=1.2, minarea=3,
                              deblend_cont=0.0001, clean='N',
                              filter_kernel=gaussian_kernel,
                              filter_type='matched',
                              segmentation_map=True)


# Visualización de fuentes detectadas
fig, ax = plt.subplots()
m, s = np.mean(weighted), np.std(weighted)
im = ax.imshow(weighted, interpolation='nearest', cmap='gray',
               vmin=m-s, vmax=m+s, origin='lower')

# plot an ellipse for each object
for i in range(len(objects)):
    e = Ellipse(xy=(objects['x'][i], objects['y'][i]),
                width=6*objects['a'][i],
                height=6*objects['b'][i],
                angle=objects['theta'][i] * 180. / np.pi)
    e.set_facecolor('none')
    e.set_edgecolor('red')
    ax.add_artist(e)

DETECTION_PARAMS = dict(
    thresh=1.2,
    minarea=3,
    kernelfwhm=3.5,
    deblend_nthresh=32,
    deblend_cont=0.0001,
    clean_param=1.0,
    clean=False,
    )

kerneldict = {}
if 'kerneltype' not in DETECTION_PARAMS.keys():
    kernel_func = Gaussian2DKernel
    kerneldict['x_stddev'] = DETECTION_PARAMS['kernelfwhm']/2.35
else:
    if DETECTION_PARAMS['kerneltype'] == 'gauss':
        kernel_func = Gaussian2DKernel
        kerneldict['x_stddev'] = DETECTION_PARAMS['kernelfwhm']/2.35

    elif DETECTION_PARAMS['kerneltype'] == 'tophat':
        kernel_func = Tophat2DKernel
        kerneldict['radius'] = DETECTION_PARAMS['kernelfwhm']


if 'kernelsize' in DETECTION_PARAMS.keys():
    kerneldict['x_size'] = DETECTION_PARAMS['kernelsize']
    kerneldict['y_size'] = DETECTION_PARAMS['kernelsize']
kerneldict['factor'] = 1
kernel = np.array(kernel_func(**kerneldict))
del DETECTION_PARAMS['kernelfwhm']
if 'kerneltype' in DETECTION_PARAMS.keys():
    del DETECTION_PARAMS['kerneltype']
if 'kernelsize' in DETECTION_PARAMS.keys():
    del DETECTION_PARAMS['kernelsize']

# Deteccion
objects1, segmap1 = sep.extract(
                weighted,
                filter_type='matched',
                filter_kernel=kernel,
                segmentation_map=True,
                **DETECTION_PARAMS
                )
fig, ax = plt.subplots()
m, s = np.mean(weighted), np.std(weighted)
im = ax.imshow(weighted, interpolation='nearest', cmap='gray',
               vmin=m-s, vmax=m+s, origin='lower')

# plot an ellipse for each object
for i in range(len(objects1)):
    e = Ellipse(xy=(objects1['x'][i], objects1['y'][i]),
                width=6*objects1['a'][i],
                height=6*objects1['b'][i],
                angle=objects1['theta'][i] * 180. / np.pi)
    e.set_facecolor('none')
    e.set_edgecolor('red')
    ax.add_artist(e)

#
w = fits.getdata('weighted_official.fits')[4303: 10045, 5280: 9455]
objects2, segmap2 = sep.extract(
                weighted,
                filter_type='matched',
                filter_kernel=kernel,
                segmentation_map=True,
                **DETECTION_PARAMS
                )
fig, ax = plt.subplots()
m, s = np.mean(w), np.std(w)
im = ax.imshow(f444w_wht, interpolation='nearest', cmap='gray',
               vmin=m-s, vmax=m+s, origin='lower')

# plot an ellipse for each object
for i in range(len(objects2)):
    e = Ellipse(xy=(objects2['x'][i], objects2['y'][i]),
                width=6*objects2['a'][i],
                height=6*objects2['b'][i],
                angle=objects2['theta'][i] * 180. / np.pi)
    e.set_facecolor('none')
    e.set_edgecolor('red')
    ax.add_artist(e)