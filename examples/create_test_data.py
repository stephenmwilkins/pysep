

import os
import sys
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import FLARE.surveys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))

import pysep.utils
import pysep.plots.image


# --- cutout a region of the HUDF to use as test data


survey_name, field_name = 'XDF', 'dXDF'

field = FLARE.surveys.surveys[survey_name].fields[field_name]

imgs = pysep.utils.images_from_field(field, verbose = True, sci_suffix = 'sci_convolved')

cutouts = {f: img.make_cutout(2000, 2000, 500) for f, img in imgs.items()}

for filter, cutout in cutouts.items():

    f = filter.split('.')[-1]

    cutout.write_to_fits(f'test_data/{f}_')


pysep.plots.image.make_flux_plots(cutouts)
