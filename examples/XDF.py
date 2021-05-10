





import os
import sys



import FLARE.surveys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))

import pysep.sep as sep
import pysep.utils
import pysep.plots.image
import pysep.analyse


survey_name, field_name = 'XDF', 'dXDF'

field = FLARE.surveys.surveys[survey_name].fields[field_name]


# --- create detection image
create_detection_image = False

if create_detection_image:
    imgs = pysep.utils.images_from_field(field, filters = field.detection_filters)
    detection_image = pysep.utils.create_stack(imgs)
    detection_image.write_to_fits(f'{field.data_dir}/detection_')



# --- read images
read_images = True

if read_images:
    # --- read individual images
    imgs = pysep.utils.images_from_field(field, verbose = True)
    # --- read detection image
    detection_image = pysep.utils.image_from_field('detection', field) # --- see "create_detection_image"


# --- run full source detection / photometry routines

run_sep = False

if run_sep: sep.run(detection_image, imgs, output_dir = field.data_dir, catalogue_name = 'cat')


# --- analyse

a = pysep.analyse.Analyser(f'{field.data_dir}/cat.h5')

a.explore_hdf5()

a.scatter('detection/xcentroid', 'detection/ycentroid', z_name = 'detection/kron_flux', logz = True)

a.scatter('detection/kron_flux', 'detection/segment_flux', logx = True, logy = True)

f = field.filters[-1]
a.scatter(f'{f}/kron_flux_nJy', f'{f}/segment_flux_nJy', logx = True, logy = True)


i = 101

x, y = a.cat['detection/xcentroid'][i], a.cat['detection/ycentroid'][i]

print('kron S/N:', a.cat['detection/kron_flux'][i]/a.cat['detection/kron_fluxerr'][i])
print('segment S/N:', a.cat['detection/segment_flux'][i]/a.cat['detection/segment_fluxerr'][i])

cutout = detection_image.make_cutout(y, x, 50)

pysep.plots.image.make_significance_plot(cutout)
