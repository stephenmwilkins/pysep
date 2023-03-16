

import os
import sys
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as cm

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))

import pysep.utils
import pysep.plots.image

# --- in this example we make an image object and plot it

data_dir = '/Users/stephenwilkins/Dropbox/Research/data/images/jwst/flags/SMACS0723/images'
filename = f'{data_dir}/jw02736-o001_t001_nircam_clear-f090w_i2d.fits'


# --- create the Image object. The Image object contains the science and weight arrays and can be easily masked
img = pysep.utils.ImageFromMultiFITS(filename)
img.measure_background_map()

# --- make a new image from a cutout of another image
cutout = img.make_cutout(500, 500, 200)


cutout.detect_sources(2.5, 11)

fig, axes = pysep.plots.image.make_mutlipanel_image(3)

axes[0] = pysep.plots.image.img_panel(axes[0], cutout.data)
axes[1] = pysep.plots.image.significance_panel(axes[1], cutout)
axes[2] = pysep.plots.image.segm_panel(axes[2], cutout.segm_detect)

plt.show()
