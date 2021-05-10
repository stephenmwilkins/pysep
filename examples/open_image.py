

import os
import sys
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import FLARE.surveys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))

import pysep.utils
import pysep.plots.image



survey_name, field_name = 'XDF', 'dXDF'

field = FLARE.surveys.surveys[survey_name].fields[field_name]

f = field.filters[-1] # final filter



img = pysep.utils.image_from_field(f, field, verbose = True)


# plt.imshow(img.sci)
# plt.show()
#
# plt.imshow(img.wht)
# plt.show()
#
# plt.imshow(img.mask)
# plt.show()


cutout = img.make_cutout(2000, 2000,500)

pysep.plots.image.make_significance_plot(cutout)
