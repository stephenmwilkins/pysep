

import os
import sys
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as cm

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))

import pysep.utils
import pysep.plots.image

# --- in this example we make an image object and plot it

filename = 'data/test/f160w'


# --- create the Image object. The Image object contains the science and weight arrays and can be easily masked
img = pysep.utils.ImageFromFITS(filename)



# --- display the science, weight, and S/N images
pysep.plots.image.make_flux_plot(img) # plot science image with linear scaling
pysep.plots.image.make_flux_plot(img, scaling = np.log10) # plot science image with log10 scaling
pysep.plots.image.make_flux_plot(img, scaling = np.arcsinh) # plot science image with arcsinh scaling

pysep.plots.image.make_flux_plot(img, ext = 'wht') # plot weight image

pysep.plots.image.make_significance_plot(img) # make nice significance image. Here the greyscale denotes pixels S/N<2 while the colour scale denotes pixels S/N>2


# --- make a new image from a cutout of another image
cutout = img.make_cutout(100, 100, 100)

pysep.plots.image.make_flux_plot(cutout) # --- plot the cutout science image # TODO: add better scaling
plt.show()
