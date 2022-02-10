
import os
import sys

import numpy as np
import photutils

import numpy as np



# --- default parameters
default_parameters = {}
default_parameters['kron_params'] = [2.5, 1]
default_parameters['threshold'] = 2.5
default_parameters['npixels'] = 5
default_parameters['nlevels'] = 32

# --- only output these columns for the individual bands. Most of the rest *should* be repitition from the detection image
default_photometry_columns = ['kron_flux','kron_fluxerr','segment_flux','segment_fluxerr']


# --- columns to exclude from the hdf5 file because they're not 1D arrays
hdf5_exclude = ['detection/sky_centroid']

def detect_sources(detection_image, parameters = default_parameters):

    segm = photutils.detect_sources(detection_image.sn(), parameters['threshold'], npixels = parameters['npixels'])
    segm_deblended = photutils.deblend_sources(detection_image.sn(), segm, npixels = parameters['npixels'], nlevels = parameters['nlevels'])
    detection_cat = photutils.SourceCatalog(detection_image.sci, segm_deblended, error = detection_image.noise, kron_params = parameters['kron_params'])

    return detection_cat, segm_deblended


def detect_sources_method(self):

    self.detection_cat, self.segm_deblended = detect_sources(self.detection_image, parameters = self.parameters)

    # --- create detetion property table
    detection_tab = self.detection_cat.to_table()

    # --- output detection property to output dictionary
    for col, colname in zip(detection_tab.itercols(), detection_tab.colnames): # extract the table into a dictionary for easier output later
        self.o[f'detection/{colname}'] = np.array(col.data)

    # --- output detection signal-to-noise
    self.o['detection/sn'] = self.o['detection/segment_flux']/self.o['detection/segment_fluxerr']

    # --- define number of sources detected
    self.N = len(self.o['detection/sn'])
    if self.verbose:
        print('-'*20)
        print(f'Number of detected sources: {self.N}')



def perform_photometry(detection_cat, segm_deblended, imgs, parameters = default_parameters):
    """ Perform the standard photometry that photutils does """
    photometry_cat = {}
    for f, img in imgs.items():
        photometry_cat[f] = photutils.SourceCatalog(img.sci, segm_deblended, error = img.noise, kron_params = parameters['kron_params'], detection_cat = detection_cat)
    return photometry_cat

def perform_photometry_method(self, photometry_columns = default_photometry_columns):

    self.photometry_cat = perform_photometry(self.detection_cat, self.segm_deblended, self.imgs, parameters = self.parameters)

    if self.verbose:
        print('-'*20)
        print('performing photometry')

    for f, photo_cat in self.photometry_cat.items():

        if self.verbose: print(f)

        if photometry_columns:
            photo_tab = photo_cat.to_table(photometry_columns) # use defined subset of columns
        else:
            photo_tab = photo_cat.to_table() # use all columns

        for col, colname in zip(photo_tab.itercols(), photo_tab.colnames):
            self.o[f'{f}/{colname}'] = np.array(col.data)

        # --- convert to nJy
        for phot_type in ['segment', 'kron']:
            self.o[f'{f}/{phot_type}_flux'] /=  self.imgs[f].nJy_to_es
            self.o[f'{f}/{phot_type}_fluxerr'] /=  self.imgs[f].nJy_to_es




def perform_aperture_photometry(detection_cat, segm_deblended, imgs, apertures = []):
    """ Perform aperture photometry -- NOT YET IMPLEMENTED """

    aperture_photometry_cat = {}
    print('WARNING: NOT YET IMPLEMENTED')

    # for f, img in imgs.items():
    #     for aperture in apertures:

    return aperture_photometry_cat

def perform_aperture_photometry_method(self):

    self.aperture_photometry_cat = perform_aperture_photometry(self.detection_cat, self.segm_deblended, self.imgs)
    print('WARNING: NOT YET IMPLEMENTED')




class SEP():

    """ The SEP class which keeps everything together and provides handy export options """

    detect_sources = detect_sources_method
    perform_photometry = perform_photometry_method
    perform_aperture_photometry = perform_aperture_photometry_method


    def __init__(self, detection_image, imgs, verbose = False, output_dir = None, parameters = default_parameters):

        self.detection_image = detection_image
        self.imgs = imgs
        self.verbose = verbose
        self.filters = list(self.imgs.keys())
        self.parameters = parameters

        if self.verbose:
            print('-'*20)
            print('list of filters:', self.filters)
            print('-'*20)
            for k,v in self.parameters.items():
                print(k, v)

        # --- output directory
        self.output_dir = output_dir
        self.o = {}


    # --- Export to Astropy tables
    def export_to_table(self, output_dir = None, cat_name = 'cat', fmt = "fits"):
        import astropy
        from astropy.io import fits, ascii
        from astropy.table import Table
        #print('Sorry not yet implemented')
        if self.verbose:
            print('-'*20)
            print('** Exporting to Astropy Table with format %s **' % fmt)
            # --- print a list of all the output columns
            for k in self.o.keys():
                print(k)
                
        if (not output_dir) & (self.output_dir is not None):
            output_dir = self.output_dir
        if (not output_dir) & (not self.output_dir):
            print('WARNING: No output directory set')
            
        t = Table([v for k,v in self.o.items() if k not in hdf5_exclude], names=[k for k,v in self.o.items() if k not in hdf5_exclude])
        print(t)
        print(t.info)
        if fmt == "ascii":
          ascii.write(t, f'{output_dir}/{cat_name}.dat', overwrite=True)
        elif fmt == "fits":
          t.write(f'{output_dir}/{cat_name}.fits', format='fits', overwrite=True)

    # --- Export to SExtractor file format.
    def export_to_SExtractor(self):
        print('Sorry not yet implemented')


    # --- Export ALL to Python pickle.
    def export_all_to_pickle(self, output_dir = None):

        """ dump everything to Python pickle. Not recommended."""

        import pickle

        pickle.dump(self, open(f'{output_dir}/full.pck','w'))

    # --- Export ALL to Python pickle.
    def export_to_pickle(self, output_dir = None, cat_name = 'cat'):

        """ dump output dictionary to Python pickle ."""

        import pickle

        pickle.dump(self, open(f'{output_dir}/{cat_name}.pck','w'))


    # --- Export to HDF5 format.
    def export_to_hdf5(self, output_dir = None, cat_name = 'cat', return_hf = False):

        import h5py

        if self.verbose:
            print('-'*20)
            print('** Exporting to HDF5 **')
            # --- print a list of all the output columns
            for k in self.o.keys():
                print(k)

        if (not output_dir) & (self.output_dir is not None):
            output_dir = self.output_dir
        if (not output_dir) & (not self.output_dir):
            print('WARNING: No output directory set')

        # --- make directory structure for the output files
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        hf = h5py.File(f'{output_dir}/{cat_name}.h5', 'w')

        # hf.attrs['N'] = self.N

        for k, v in self.o.items():
            print(k)
            if k not in hdf5_exclude:
                hf.create_dataset(k, data = v)

        if return_hf:
            return hf
        else:
            hf.flush()


    def run_default(self):

        self.detect_sources()
        self.perform_photometry()
        self.perform_aperture_photometry()
        self.export_to_hdf5()
