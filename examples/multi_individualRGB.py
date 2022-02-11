

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

filename = lambda f: f'{data_dir}/{f}' # define a function that gives the full filename for different filters

imgs = {f: pysep.utils.ImageFromFITS(filename(f)) for f in filters}

r = pysep.utils.create_stack({f:imgs[f] for f in ['f105w','f125w','f140w','f160w']})
g = pysep.utils.create_stack({f:imgs[f] for f in ['f775w','f814w','f850lp']})
b = pysep.utils.create_stack({f:imgs[f] for f in ['f435w','f606w']})

r.data[r.data!=r.data] = 0.
g.data[g.data!=g.data] = 0.
b.data[b.data!=b.data] = 0.

# r.data /= np.sum(r.data)
# g.data /= np.sum(g.data)
# b.data /= np.sum(b.data)

r.data /= np.max(r.data)
g.data /= np.max(g.data)
b.data /= np.max(b.data)


detection_image = imgs['f160w'] # simply use the f160w band

# --- initialise SEP, this keeps everything together and provides handy output methods

SEP = sep.SEP(detection_image, imgs, verbose = True) # you could at this point change the default parameters

SEP.detect_sources() # detect sources


# for i in range(SEP.N):
for i in range(5):

    x, y = SEP.o['detection/xcentroid'][i], SEP.o['detection/ycentroid'][i]

    # --- make RGB image of the galaxy

    r_, g_, b_ = [img.make_cutout(y, x, 100, calculate_background = False).data for img in [r,g,b]] # make cutouts in each bands

    rgb_ = r_ + g_ + b_

    r_ /= np.max(rgb_)/3
    g_ /= np.max(rgb_)/3
    b_ /= np.max(rgb_)/3


    pysep.plots.image.make_rgb_image(r_,g_,b_)
