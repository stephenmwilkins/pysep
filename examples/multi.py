

import os
import sys

import numpy as np

import pysep.sep as sep
import pysep.utils
import pysep.plots.image
import pysep.analyse



# --- read the test data. The test data is organised as: data/test/{f}_{image_type}.fits where f is the filter and image_type is e.g. sci, wht


# --- default parameters
parameters = {}
parameters['kron_params'] = [2.5, 1]
parameters['npixels'] = 5
parameters['nlevels'] = 32
parameters['nsigma'] = 2
parameters['deblend_contrast'] = 0.001

parameters['smooth'] = False
# parameters['smooth'] = {'smooth_fwhm':2, 'kernel_size':5}



data_dir = 'data/test'

filters = ['f435w','f606w','f775w','f814w','f850lp','f105w','f125w','f140w','f160w']

zeropoints = {'f105w': 26.269, 'f125w': 26.230, 'f140w': 26.452, 'f160w': 25.946, 'f435w': 25.684, 'f606w': 26.505, 'f775w': 25.678, 'f814w': 25.959, 'f850lp': 24.867}


filename = lambda f, ext: f'{data_dir}/{f}_{ext}' # define a function that gives the full filename for different filters

# --- create a dictionary of images
imgs = {f: pysep.utils.ImageFromFITS(filename(f,'sci'), filename(f,'wht'), zeropoint = zeropoints[f]) for f in filters}

# pysep.plots.image.make_significance_plots(imgs) # make a nice S/N plot of all filters


# --- for PYSEP to work we need to define (or create) a detection image

detection_image = imgs['f160w'] # simply use the f160w band
# detection_image = pysep.utils.ImageFromFITS(filename('detection')) # use a separate detection image in the same folder

# alternatively we can create a detection image from a weighted stack of a set of images
# detection_filters = ['f105w','f125w','f140w','f160w'] # all WFC3 images
# detection_image = pysep.utils.create_stack({f:imgs[f] for f in detection_filters})

pysep.plots.image.make_significance_plot(detection_image) # plot significance map



# --- initialise SEP, this keeps everything together and provides handy output methods

SEP = sep.SEP(verbose = True, parameters = parameters) # you could at this point change the default parameters

SEP.detect_sources(detection_image) # detect sources

# pysep.plots.image.make_segm_plot(SEP.segm_deblended) # plot segmentation map

SEP.perform_photometry(imgs) # perform matched kron and segment photometry on all images using positions and apertures based on the detection images


# SEP.export_to_pickle(output_dir = data_dir) # dump output dictionary to Python pickle
# SEP.export_to_hdf5(output_dir = data_dir) # dump output dictionary to HDF5
# SEP.export_to_pandas(output_dir = data_dir) # dump output dictionary to Python pickle
# SEP.export_to_table(output_dir = data_dir) # dump output dictionary to astropy table. NOT IMPLEMENTED
# SEP.export_to_SExtractor(output_dir = data_dir) # dump output dictionary to SExtractor format. NOT IMPLEMENTED
