


import numpy as np

from astropy.io import fits

from photutils import Background2D, MedianBackground

# from photutils import CircularAperture
# from photutils import aperture_photometry





class empty: pass

def make_cutout(data, x, y, width):

    """extract cut out from arbitrary data"""

    cutout = np.zeros((width, width))

    xmin = x - width // 2
    xmax = x + width // 2
    ymin = y - width // 2
    ymax = y + width // 2

    xstart = 0
    ystart = 0
    xend = width
    yend = width

    if xmin < 0:
        xstart = -xmin
        xmin = 0
    if ymin < 0:
        ystart = -ymin
        ymin = 0
    if xmax > data.shape[0]:
        xend -= xmax - data.shape[0]
        xmax = data.shape[0]
    if ymax > data.shape[1]:
        yend -= ymax - data.shape[1]
        ymax = data.shape[1]

    data = np.array(data)
    cutout[xstart:xend,ystart:yend] = data[xmin:xmax,ymin:ymax]

    return cutout




# --- this is the Image class which is extensively used

class Image:

    def measure_background_and_rms(self):

        self.measure_background_map()

        self.background_rms = 1 / np.sqrt(self.wht)

        if self.exposure_time:

            self.exposure_time_map = self.exposure_time * self.background_map.background_rms_median**2 * self.wht

            # effective gain parameter required to be positive everywhere (not zero), so adding small value 1e-8
            self.data_rms = calc_total_error(self.data, self.background_rms, self.exposure_time_map+1e-8)

        else:

            print('WARNING: Using naive data_rms calculation')
            self.data_rms = self.background_rms



    def sn(self):

        return self.data/self.data_rms


    def get_random_location(self):

        """get (single) random location on the image"""

        pos = np.random.choice(self.data.count())
        return np.take((~self.data.mask).nonzero(), pos, axis=1)

    def get_random_locations(self, N):

        """get N random locations on the image"""

        pos = np.random.choice(self.data.count(), size=N)
        return np.take((~self.data.mask).nonzero(), pos, axis=1)


    def get_area(self):

        """calculated non-masked area in units of arcmin2"""

        return self.data.count()*self.pixel_scale**2/3600.


    def make_cutout(self, x, y, width, calculate_background = True):

        """extract cut out"""

        data = np.zeros((width, width))
        wht = np.zeros((width, width))

        x = int(np.round(x, 0))
        y = int(np.round(y, 0))

        xmin = x - width // 2
        xmax = x + width // 2
        ymin = y - width // 2
        ymax = y + width // 2

        xstart = 0
        ystart = 0
        xend = width
        yend = width

        if xmin < 0:
            xstart = -xmin
            xmin = 0
        if ymin < 0:
            ystart = -ymin
            ymin = 0
        if xmax > self.data.shape[0]:
            xend -= xmax - self.data.shape[0]
            xmax = self.data.shape[0]
        if ymax > self.data.shape[1]:
            yend -= ymax - self.data.shape[1]
            ymax = self.data.shape[1]

#         if self.verbose: print(xmin, xmax, ymin, ymax)
#         if self.verbose: print(xstart, xend)
#         if self.verbose: print(ystart, yend)
#         if self.verbose: print(data.shape)
#         if self.verbose: print(self.data[xmin:xmax,ymin:ymax].shape)

        if (width % 2) != 0:
            xmax += 1
            ymax += 1

        data[xstart:xend,ystart:yend] = self.data[xmin:xmax,ymin:ymax]
        wht[xstart:xend,ystart:yend] = self.wht[xmin:xmax,ymin:ymax]

        return ImageFromArrays(data, wht, self.pixel_scale, zeropoint = self.zeropoint, nJy_to_es = self.nJy_to_es, verbose = self.verbose, calculate_background = calculate_background)



    def determine_depth(self, N = 10000, aperture_diameter_arcsec = 0.35, sigma = 5.):

        """determine depth using random apertures"""

        aperture_centres = tuple(self.get_random_locations(N).T)
        apertures = [CircularAperture(aperture_centres, r=r) for r in [(aperture_diameter_arcsec/self.pixel_scale)/2.]] # r in pixels
        phot_table = aperture_photometry(self.data, apertures)
        aperture_fluxes = phot_table['aperture_sum_0'].quantity
        negative_aperture_fluxes = aperture_fluxes[aperture_fluxes<0]
        return -np.percentile(negative_aperture_fluxes, 100.-68.3) * sigma


    def write_to_fits(self, filename = 'temp/'):

        sci_hdu = fits.PrimaryHDU(self.data)
        sci_hdu.writeto(f'{filename}sci.fits')

        wht_hdu = fits.PrimaryHDU(self.wht)
        wht_hdu.writeto(f'{filename}wht.fits')

        rms_hdu = fits.PrimaryHDU(self.data_rms)
        rms_hdu.writeto(f'{filename}data_rms.fits')


    def measure_background_map(self, bkg_size=50, filter_size=3, verbose=True):
        # Calculate sigma-clipped background in cells of 50x50 pixels, then median filter over 3x3 cells
        # For best results, the image should span an integer number of cells in both dimensions (here, 1000=20x50 pixels)
        # https://photutils.readthedocs.io/en/stable/background.html
        self.background_map = Background2D(self.data, bkg_size, filter_size = filter_size)



# --- designed to work with Steve's FLARE module. Omitted here.
# def images_from_field(field, filters = None, verbose = False, sci_suffix = 'sci', wht_suffix = 'wht'):
#
#     # --- uses a flare.surveys field object to set the relevant parameters
#
#     if field.mask_file:
#         mask = fits.getdata(f'{field.data_dir}/{field.mask_file}')
#     else:
#         mask = None
#
#     if not filters:
#         filters = field.filters
#
#     return {filter: ImageFromFile(field.data_dir, filter, mask = mask, pixel_scale = field.pixel_scale, verbose = verbose, sci_suffix = sci_suffix, wht_suffix = wht_suffix) for filter in filters}



# --- designed to work with Steve's FLARE module. Omitted here.
# def image_from_field(filter, field, verbose = False, sci_suffix = 'sci', wht_suffix = 'wht'):
#
#     # --- uses a flare.surveys field object to set the relevant parameters
#
#     if field.mask_file:
#         mask = fits.getdata(f'{field.data_dir}/{field.mask_file}')
#     else:
#         mask = None
#
#     return ImageFromFile(field.data_dir, filter, mask = mask, pixel_scale = field.pixel_scale, verbose = verbose, sci_suffix = sci_suffix, wht_suffix = wht_suffix)



# --- this class is used to read in an image from a FITS file

class ImageFromFITS(Image):

    def __init__(self, filename, filter = None, mask = None, pixel_scale = 0.06, verbose = False, sci_suffix = 'sci', wht_suffix = 'wht', zeropoint = None, nJy_to_es = None, exposure_time = None):

        """generate instance of image class from file"""

        self.verbose = verbose

        self.exposure_time = exposure_time

        self.filter = filter
        self.pixel_scale = pixel_scale

        # self.sci = fits.getdata(f'{data_dir}/{f}_{sci_suffix}.fits')
        # self.wht = fits.getdata(f'{data_dir}/{f}_{wht_suffix}.fits')

        self.data = fits.getdata(f'{filename}_{sci_suffix}.fits')
        self.wht = fits.getdata(f'{filename}_{wht_suffix}.fits')


        # --- define image zeropoint and/or conversion from nJy to electrons per second

        if (zeropoint is None) and (nJy_to_es is None): print('WARNING: no zeropoint set. This is needed for photometry.')

        if zeropoint and (nJy_to_es is None):
            nJy_to_es = 1E-9 * 10**(0.4*(zeropoint-8.9))

        self.zeropoint = zeropoint # AB magnitude zeropoint
        self.nJy_to_es = nJy_to_es # conversion from nJy to e/s

        self.mask = mask

        if type(mask) == np.ndarray:
            self.mask = mask
        else:
            self.mask = (self.wht == 0)

        self.data = np.ma.masked_array(self.data, mask = self.mask)
        self.wht = np.ma.masked_array(self.wht, mask = self.mask)

        if verbose:
            print(f'shape: ', self.data.shape)

        self.measure_background_and_rms()



# --- this class is used to read in an image from a FITS file


class ImageFromArrays(Image):

    def __init__(self, data, wht, pixel_scale, zeropoint = False, nJy_to_es = False,  verbose = False, exposure_time = None, calculate_background = True):

        """generate instance of image class from cutout"""

        self.verbose = verbose

        self.exposure_time = exposure_time

        self.pixel_scale = pixel_scale
        self.zeropoint = zeropoint # AB magnitude zeropoint
        self.nJy_to_es = nJy_to_es # conversion from nJy to e/s

        self.data = data
        self.wht = wht

        if calculate_background:
            self.measure_background_and_rms()





def create_stack(imgs):

    first_img = next(iter(imgs.values()))

    shape = first_img.data.shape
    data = np.zeros(shape)
    wht = np.zeros(shape)

    for filter, img in imgs.items():
        data += img.data * img.wht
        wht += img.wht

    data /= wht

    return ImageFromArrays(data, wht, first_img.pixel_scale)
