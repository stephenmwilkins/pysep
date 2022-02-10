

import os
import sys
import numpy as np
import photutils

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
import pysep.sep as sep
import pysep.utils
# import pysep.analyse
import pysep.plots.image


# --- this is verbose exploration of what PYSEP will automatically. PYSEP will however do this for multiple filters/bands!


filename = 'data/test/f160w'

# --- create the Image object. The Image object contains the science and weight arrays and can be easily masked
img = pysep.utils.ImageFromFITS(filename)

# --- make nice S/N image. Here the greyscale denotes pixels S/N<2 while the colour scale denotes pixels S/N>2
pysep.plots.image.make_significance_plot(img)

# --- create the segmention image
detection_image = img.sn() # define a signal-to-noise image which will be used as the detection image
threshold = 2.5 # the S/N threshold needed for pixels
npixels = 5 # the minium number of connected pixels
segm = photutils.detect_sources(detection_image, threshold, npixels = npixels) # create segmentation image
pysep.plots.image.make_segm_plot(segm) # make nice segmentation image plot

# --- create the DE-BLENDED segmention image
nlevels = 32 # number of deblending levels
segm_deblended = photutils.deblend_sources(detection_image, segm, npixels = npixels, nlevels = nlevels)
pysep.plots.image.make_segm_plot(segm) # make nice deblended segmentation image plot

# --- create a catalog of photometry and morphological properties for sources defined by a segmentation image.
kron_params = [2.5, 1]
detection_cat = photutils.SourceCatalog(img.sci, segm_deblended, error = img.noise, kron_params = kron_params)
detection_tab = detection_cat.to_table() # convert the detection_cat to a detection)table
for colname in detection_tab.colnames: print(colname) # print column names to see what we have

# --- perform aperture photometry on the above (as it's not done by default)
