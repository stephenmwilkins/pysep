


import numpy as np

from astropy.io import fits
# from photutils import CircularAperture
# from photutils import aperture_photometry

#import flare
#import flare.observatories

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


    def sn(self):

        return self.sci/self.noise


    def get_random_location(self):

        """get (single) random location on the image"""

        pos = np.random.choice(self.sci.count())
        return np.take((~self.sci.mask).nonzero(), pos, axis=1)

    def get_random_locations(self, N):

        """get N random locations on the image"""

        pos = np.random.choice(self.sci.count(), size=N)
        return np.take((~self.sci.mask).nonzero(), pos, axis=1)


    def get_area(self):

        """calculated non-masked area in units of arcmin2"""

        return self.sci.count()*self.pixel_scale**2/3600.


    def make_cutout(self, x, y, width):

        """extract cut out"""

        sci = np.zeros((width, width))
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
        if xmax > self.sci.shape[0]:
            xend -= xmax - self.sci.shape[0]
            xmax = self.sci.shape[0]
        if ymax > self.sci.shape[1]:
            yend -= ymax - self.sci.shape[1]
            ymax = self.sci.shape[1]

#         if self.verbose: print(xmin, xmax, ymin, ymax)
#         if self.verbose: print(xstart, xend)
#         if self.verbose: print(ystart, yend)
#         if self.verbose: print(sci.shape)
#         if self.verbose: print(self.sci[xmin:xmax,ymin:ymax].shape)

        if (width % 2) != 0:
            xmax += 1
            ymax += 1

        sci[xstart:xend,ystart:yend] = self.sci[xmin:xmax,ymin:ymax]
        wht[xstart:xend,ystart:yend] = self.wht[xmin:xmax,ymin:ymax]

        return ImageFromArrays(sci, wht, self.pixel_scale, zeropoint = self.zeropoint, nJy_to_es = self.nJy_to_es, verbose = self.verbose)



    def determine_depth(self, N = 10000, aperture_diameter_arcsec = 0.35, sigma = 5.):

        """determine depth using random apertures"""

        aperture_centres = tuple(self.get_random_locations(N).T)
        apertures = [CircularAperture(aperture_centres, r=r) for r in [(aperture_diameter_arcsec/self.pixel_scale)/2.]] # r in pixels
        phot_table = aperture_photometry(self.sci, apertures)
        aperture_fluxes = phot_table['aperture_sum_0'].quantity
        negative_aperture_fluxes = aperture_fluxes[aperture_fluxes<0]
        return -np.percentile(negative_aperture_fluxes, 100.-68.3) * sigma


    def write_to_fits(self, filename = 'temp/'):

        sci_hdu = fits.PrimaryHDU(self.sci)
        sci_hdu.writeto(f'{filename}sci.fits')

        wht_hdu = fits.PrimaryHDU(self.wht)
        wht_hdu.writeto(f'{filename}wht.fits')

        rms_hdu = fits.PrimaryHDU(self.noise)
        rms_hdu.writeto(f'{filename}rms.fits')





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

    def __init__(self, filename, filter = None, mask = None, pixel_scale = 0.06, verbose = False, sci_suffix = 'sci', wht_suffix = 'wht', zeropoint = None, nJy_to_es = None):

        """generate instance of image class from file"""

        self.verbose = verbose

        self.filter = filter
        self.pixel_scale = pixel_scale

        # self.sci = fits.getdata(f'{data_dir}/{f}_{sci_suffix}.fits')
        # self.wht = fits.getdata(f'{data_dir}/{f}_{wht_suffix}.fits')

        self.sci = fits.getdata(f'{filename}_{sci_suffix}.fits')
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

        self.sci = np.ma.masked_array(self.sci, mask = self.mask)
        self.wht = np.ma.masked_array(self.wht, mask = self.mask)

        if verbose:
            print(f'shape: ', self.sci.shape)

        self.noise = 1./np.sqrt(self.wht)
        # self.sig = self.sci/self.noise



# --- this class is used to read in an image from a FITS file


class ImageFromArrays(Image):

    def __init__(self, sci, wht, pixel_scale, zeropoint = False, nJy_to_es = False,  verbose = False):

        """generate instance of image class from cutout"""

        self.verbose = verbose

        self.pixel_scale = pixel_scale
        self.zeropoint = zeropoint # AB magnitude zeropoint
        self.nJy_to_es = nJy_to_es # conversion from nJy to e/s

        self.sci = sci
        self.wht = wht
        self.noise = 1./np.sqrt(self.wht)
        # self.sig = self.sci/self.noise



def create_stack(imgs):

    first_img = next(iter(imgs.values()))

    shape = first_img.sci.shape
    sci = np.zeros(shape)
    wht = np.zeros(shape)

    for filter, img in imgs.items():
        sci += img.sci * img.wht
        wht += img.wht

    sci /= wht

    return ImageFromArrays(sci, wht, first_img.pixel_scale)
