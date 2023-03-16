


import astropy  # version 4.2 is required to write magnitudes to ecsv file
from astropy.io import fits

data_dir = '/Users/stephenwilkins/Dropbox/Research/data/images/jwst/flags/SMACS0723/images'
filename = f'{data_dir}/jw02736-o001_t001_nircam_clear-f090w_i2d.fits'

data_dir = '/Users/stephenwilkins/Dropbox/Research/data/images/jwst/ceers/images'
filename = f'{data_dir}/ceers_nircam1_f200w_v0.2_i2d.fits'


hdul = fits.open(filename)

hdul.info()
for i,k in hdul[1].header.items():
    print(i,k)




# imwcs = wcs.WCS(hdu[0].header, hdu)
#
# weight = fits.open(weight_files[filt])[0].data
