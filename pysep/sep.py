

import numpy as np
import photutils
import h5py
import numpy as np

import FLARE.observatories




kron_params = [2.5, 1]
threshold = 2.5
npixels = 5
nlevels = 32
photometry_columns = ['kron_flux', 'kron_fluxerr', 'segment_flux', 'segment_fluxerr']
detection_columns = photometry_columns + ['kron_radius', 'xcentroid', 'ycentroid', 'area', 'bbox_xmax', 'bbox_xmin', 'bbox_ymax', 'bbox_ymin', 'ellipticity', 'equivalent_radius', 'gini']


def detect_sources(detection_image):

    segm = photutils.detect_sources(detection_image.sn(), threshold, npixels = npixels)
    segm_deblended = photutils.deblend_sources(detection_image.sn(), segm, npixels = npixels, nlevels = nlevels)
    detection_cat = photutils.SourceCatalog(detection_image.sci, segm_deblended, error = detection_image.noise, kron_params = kron_params)

    print(len(detection_cat))

    return detection_cat, segm_deblended


def perform_photometry(detection_cat, segm_deblended, imgs):

    photometry_cat = {}
    for f, img in imgs.items():
        photometry_cat[f] = photutils.SourceCatalog(img.sci, segm_deblended, error = img.noise, kron_params = kron_params, detection_cat = detection_cat)

    return photometry_cat


def save_as_hdf5(output_dir, detection_cat, photometry_cat = None, catalogue_name = 'cat'):

    hf = h5py.File(f'{output_dir}/{catalogue_name}.h5', 'w')

    detection_tab = detection_cat.to_table(detection_columns)

    for col, colname in zip(detection_tab.itercols(), detection_tab.colnames):

        hf[f'detection/{colname}'] = np.array(col.data)


    if photometry_cat:

        for f, photo_cat in photometry_cat.items():

            photo_tab = photo_cat.to_table(photometry_columns)

            for col, colname in zip(photo_tab.itercols(), photo_tab.colnames):

                hf[f'{f}/{colname}'] = np.array(col.data)

            for phot_type in ['segment', 'kron']:

                hf[f'{f}/{phot_type}_flux_nJy'] = hf[f'{f}/{phot_type}_flux'][:] / FLARE.observatories.filter_info[f]['nJy_to_es']
                hf[f'{f}/{phot_type}_fluxerr_nJy'] = hf[f'{f}/{phot_type}_fluxerr'][:] / FLARE.observatories.filter_info[f]['nJy_to_es']



    hf.close()



def run(detection_image, imgs, output_dir = None, catalogue_name = 'cat'):

    detection_cat, segm_deblended = detect_sources(detection_image)

    photometry_cat = perform_photometry(detection_cat, segm_deblended, imgs)

    if output_dir:

        save_as_hdf5(output_dir, detection_cat, photometry_cat = photometry_cat, catalogue_name = catalogue_name)


class Runner():

    # --- keeps hold of everything

    def __init__(self, detection_image, imgs, output_dir = None, catalogue_name = 'cat'):

        self.detection_image = detection_image
        self.imgs = imgs
        self.detection_cat, self.segm_deblended = detect_sources(self.detection_image)
        self.photometry_cat = perform_photometry(self.detection_cat, self.segm_deblended, self.imgs)

        if output_dir:
            save_as_hdf5(output_dir, self.detection_cat, photometry_cat = self.photometry_cat, catalogue_name = catalogue_name)
