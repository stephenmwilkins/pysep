

import os
import sys
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as cm

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))

import pysep.utils
import pysep.plots.image

# --- in this example we make an image object and plot it

# data_dir = '/Users/stephenwilkins/Dropbox/Research/data/images/jwst/flags/SMACS0723/images'
# filename = f'{data_dir}/jw02736-o001_t001_nircam_clear-f090w_i2d.fits'

data_dir = '/Users/stephenwilkins/Dropbox/Research/data/images/jwst/ceers/images'
filename = f'{data_dir}/ceers_nircam1_f200w_v0.2_i2d.fits'

# --- create the Image object. The Image object contains the science and weight arrays and can be easily masked
img = pysep.utils.ImageFromMultiFITS(filename)
img.measure_background_map()

# --- display the science, weight, and S/N images
# pysep.plots.image.make_flux_plot(img) # plot science image with linear scaling
# pysep.plots.image.make_flux_plot(img, scaling = np.log10) # plot science image with log10 scaling
# pysep.plots.image.make_flux_plot(img, scaling = np.arcsinh) # plot science image with arcsinh scaling
# pysep.plots.image.make_flux_plot(img, ext = 'wht') # plot weight image

# pysep.plots.image.make_significance_plot(img) # make nice significance image. Here the greyscale denotes pixels S/N<2 while the colour scale denotes pixels S/N>2


for im in [img.data, img.wht, img.bkg, img.bkg_rms]:
    print(im.shape, np.median(im), np.mean(im), np.std(im), im.dtype)


# --- make a new image from a cutout of another image
cutout = img.make_cutout(500, 500, 200)


fig, ax = pysep.plots.image.make_image_plot(cutout.data) # --- plot the cutout science image # TODO: add better scaling
plt.show()

# fig, ax = pysep.plots.image.make_image_plot(cutout.data, scaling = np.log10) # --- plot the cutout science image # TODO: add better scaling
# plt.show()
#
# fig, ax = pysep.plots.image.make_image_plot(cutout.data, scaling = np.arcsinh) # --- plot the cutout science image # TODO: add better scaling
# plt.show()

fig, ax = pysep.plots.image.make_images_plot([cutout.data, cutout.err, cutout.wht, cutout.bkg, cutout.bkg_rms]) # --- plot the cutout science image # TODO: add better scaling
plt.show()

fig, ax = pysep.plots.image.make_significance_plot(cutout)
plt.show()
