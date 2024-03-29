

import os
import sys

import numpy as np

import pysep.sep as sep
import pysep.utils
import pysep.plots.image
import pysep.analyse



# --- read the test data. The test data is organised as: data/test/{f}_{image_type}.fits where f is the filter and image_type is e.g. sci, wht


data_dir = 'data/test'

filters = ['f435w','f606w','f775w','f814w','f850lp','f105w','f125w','f140w','f160w']

filename = lambda f, ext: f'{data_dir}/{f}_{ext}' # define a function that gives the full filename for different filters

# --- create a dictionary of images
imgs = {f: pysep.utils.ImageFromFITS(filename(f,'sci'), filename(f,'wht'), zeropoint = zeropoints[f]) for f in filters}

# --- for PYSEP to work we need to define (or create) a detection image

detection_image = imgs['f160w'] # simply use the f160w band
# detection_image = pysep.utils.ImageFromFITS(filename('detection')) # use a separate detection image in the same folder

# alternatively we can create a detection image from a weighted stack of a set of images
# detection_filters = ['f105w','f125w','f140w','f160w'] # all WFC3 images
# detection_image = pysep.utils.create_stack({f:imgs[f] for f in detection_filters})

# --- initialise SEP, this keeps everything together and provides handy output methods

SEP = sep.SEP(verbose = True) # you could at this point change the default parameters

SEP.detect_sources(detection_image) # detect sources

SEP.perform_photometry(imgs) # perform matched kron and segment photometry on all images using positions and apertures based on the detection images



# for i in range(SEP.N):
for i in range(4):

    print('-'*20)
    print('Exploring individual galaxies')

    i = 4 # choose a random galaxy

    x, y = SEP.o['detection/xcentroid'][i], SEP.o['detection/ycentroid'][i]

    print(f"kron S/N: {SEP.o['detection/kron_flux'][i]/SEP.o['detection/kron_fluxerr'][i]:.1f}")
    print(f"segment S/N: {SEP.o['detection/segment_flux'][i]/SEP.o['detection/segment_fluxerr'][i]:.1f}")

    # --- make a cutout of the galaxy and place the Kron aperture on top
    cutout = detection_image.make_cutout(y, x, 50)
    ap = SEP.detection_cat.kron_aperture[i]
    import matplotlib.patches as mpatches
    theta_deg = ap.theta * 180. / np.pi
    patch = mpatches.Ellipse([25,25], 2.*ap.a, 2.*ap.b, theta_deg)
    pysep.plots.image.make_significance_plot(cutout, patches = [patch])
