

import os
import sys

import numpy as np

import pysep.sep as sep
import pysep.utils
import pysep.plots.image
import pysep.analyse



# --- read the test data. The test data is organised as: data/test/{f}_{image_type}.fits where f is the filter and image_type is e.g. sci, wht


data_dir = 'examples/data/test'

filters = ['f435w','f606w','f775w','f814w','f850lp','f105w','f125w','f140w','f160w']

zeropoints = {'f105w': 26.269, 'f125w': 26.230, 'f140w': 26.452, 'f160w': 25.946, 'f435w': 25.684, 'f606w': 26.505, 'f775w': 25.678, 'f814w': 25.959, 'f850lp': 24.867}


filename = lambda f: f'{data_dir}/{f}' # define a function that gives the full filename for different filters

# --- create a dictionary of images
imgs = {f: pysep.utils.ImageFromFITS(filename(f), zeropoint = zeropoints[f]) for f in filters}

pysep.plots.image.make_significance_plots(imgs) # make a nice S/N plot of all filters



# --- for PYSEP to work we need to define (or create) a detection image

# detection_image = imgs['f160w'] # simply use the f160w band
# detection_image = pysep.utils.ImageFromFITS(filename('detection')) # use a separate detection image in the same folder

# alternatively we can create a detection image from a weighted stack of a set of images
detection_filters = ['f105w','f125w','f140w','f160w'] # all WFC3 images
detection_image = pysep.utils.create_stack({f:imgs[f] for f in detection_filters})

pysep.plots.image.make_significance_plot(detection_image) # plot significance map



# --- initialise SEP, this keeps everything together and provides handy output methods

SEP = sep.SEP(detection_image, imgs, verbose = True) # you could at this point change the default parameters

SEP.detect_sources() # detect sources

pysep.plots.image.make_segm_plot(SEP.segm_deblended) # plot segmentation map

SEP.perform_photometry() # perform matched kron and segment photometry on all images using positions and apertures based on the detection images


# SEP.export_to_pickle(output_dir = data_dir) # dump output dictionary to Python pickle
#SEP.export_to_hdf5(output_dir = data_dir) # dump output dictionary to HDF5
# SEP.export_to_pandas(output_dir = data_dir) # dump output dictionary to Python pickle
SEP.export_to_table(output_dir = data_dir) # dump output dictionary to astropy table. NOT IMPLEMENTED
# SEP.export_to_SExtractor(output_dir = data_dir) # dump output dictionary to SExtractor format. NOT IMPLEMENTED


# see analyse.py for analysing this output

# before we end let's have a look at some individual objects

explore_individual_objects = True
if explore_individual_objects:

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
