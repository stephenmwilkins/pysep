

import os
import sys

import numpy as np

import FLARE.surveys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))

import pysep.sep as sep
import pysep.utils
import pysep.plots.image
import pysep.analyse


survey_name, field_name = 'test', 'test'

field = FLARE.surveys.surveys[survey_name].fields[field_name]

imgs = pysep.utils.images_from_field(field, verbose = True)

detection_image = pysep.utils.image_from_field('detection', field) # --- see "create_detection_image"


# pysep.plots.image.make_flux_plots(imgs) # --- need to scale differently
# pysep.plots.image.make_significance_plots(imgs)
# pysep.plots.image.make_significance_plot(detection_image)


# pysep.plots.image.make_significance_segm_plot(detection_image, segm_deblended)


# sep.run(detection_image, imgs, output_dir = field.data_dir, catalogue_name = 'cat')
r = sep.Runner(detection_image, imgs, output_dir = field.data_dir, catalogue_name = 'cat') # class based alternative to the above, allows access to intermediate data




a = pysep.analyse.Analyser(f'{field.data_dir}/cat.h5')

a.explore_hdf5()

# a.scatter('detection/xcentroid', 'detection/ycentroid', z_name = 'detection/kron_flux', logz = True)
#
# a.scatter('detection/kron_flux', 'detection/segment_flux', logx = True, logy = True)
#
#
# f = field.filters[-1]
# a.scatter(f'{f}/kron_flux_nJy', f'{f}/segment_flux_nJy', logx = True, logy = True)


i = 4

x, y = a.cat['detection/xcentroid'][i], a.cat['detection/ycentroid'][i]

print('kron S/N:', a.cat['detection/kron_flux'][i]/a.cat['detection/kron_fluxerr'][i])
print('segment S/N:', a.cat['detection/segment_flux'][i]/a.cat['detection/segment_fluxerr'][i])

cutout = detection_image.make_cutout(y, x, 50)

pysep.plots.image.make_significance_plot(cutout)


ap = r.detection_cat.kron_aperture[i]

import matplotlib.patches as mpatches

theta_deg = ap.theta * 180. / np.pi

patch = mpatches.Ellipse([25,25], 2.*ap.a, 2.*ap.b, theta_deg)



pysep.plots.image.make_significance_plot(cutout, patches = [patch])
