
import os
import sys

import FLARE.surveys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))

import pysep.utils
import pysep.plots.image



survey_name, field_name = 'test', 'test'

field = FLARE.surveys.surveys[survey_name].fields[field_name]

print(field.detection_filters)

imgs = pysep.utils.images_from_field(field, filters = field.detection_filters)

pysep.plots.image.make_significance_plots(imgs)

detection_image = pysep.utils.create_stack(imgs)

pysep.plots.image.make_significance_plot(detection_image)

detection_image.write_to_fits('test_data/detection_')
