


import numpy as np

from astropy.io import fits
import astropy.wcs as wcs
import astropy.units as u
from astropy.visualization import make_lupton_rgb, SqrtStretch, LogStretch, hist, simple_norm
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.coordinates import SkyCoord


from photutils import Background2D, MedianBackground, detect_sources, deblend_sources
from photutils.utils import calc_total_error

# from photutils import CircularAperture
# from photutils import aperture_photometry

from photutils.psf.matching import resize_psf, create_matching_kernel, CosineBellWindow
from astropy.convolution import convolve, convolve_fft # , Gaussian2DKernel, Tophat2DKernel

from scipy import ndimage

from inspect import signature


JWST_flux_units = u.MJy / u.sr

output_properties = 'label xcentroid ycentroid sky_centroid area semimajor_sigma semiminor_sigma'.split()
output_properties += 'fwhm ellipticity orientation gini'.split()
output_properties += 'kron_radius local_background segment_flux segment_fluxerr kron_flux kron_fluxerr'.split()


def fluxes2mags(flux, fluxerr):
    nondet = flux < 0  # Non-detection if flux is negative
    unobs = (fluxerr <= 0) + (fluxerr == np.inf)  # Unobserved if flux uncertainty is negative or infinity

    mag = flux.to(u.ABmag)
    magupperlimit = fluxerr.to(u.ABmag) # 1-sigma upper limit if flux=0

    mag = np.where(nondet, 99 * u.ABmag, mag)
    mag = np.where(unobs, -99 * u.ABmag, mag)

    magerr = 2.5 * np.log10(1 + fluxerr/flux)
    magerr = magerr.value * u.ABmag

    magerr = np.where(nondet, magupperlimit, magerr)
    magerr = np.where(unobs, 0 * u.ABmag, magerr)

    return mag, magerr


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

    def measure_background_map(self, bkg_size=50, filter_size=3, verbose=True):
        # Calculate sigma-clipped background in cells of 50x50 pixels, then median filter over 3x3 cells
        # For best results, the image should span an integer number of cells in both dimensions (e.g., 1000=20x50 pixels)
        # https://photutils.readthedocs.io/en/stable/background.html
        self.background_map = Background2D(self.data, bkg_size, filter_size=filter_size)
        self.bkg = self.background_map.background
        self.bkg_rms = self.background_map.background_rms

    def smooth_data(self, smooth_fwhm=2, kernel_size=5):
        # convolve data with Gaussian
        # convolved_data used for source detection and to calculate source centroids and morphological properties
        smooth_sigma = smooth_fwhm * gaussian_fwhm_to_sigma
        self.smooth_kernel = Gaussian2DKernel(smooth_sigma, x_size=kernel_size, y_size=kernel_size)
        self.smooth_kernel.normalize()
        self.convolved_data = convolve(self.data, self.smooth_kernel)


    def detect_sources(self, nsigma, npixels, smooth_fwhm=2, kernel_size=5, deblend_levels=32, deblend_contrast=0.005):

        # Set detection threshold map as nsigma times RMS above background pedestal
        detection_threshold = (nsigma * self.bkg_rms) + self.bkg
        # detection_threshold = (nsigma * self.background_map.background_rms)/self.wht + self.background_map.background

        # Before detection, convolve data with Gaussian
        self.smooth_data(smooth_fwhm, kernel_size)

        # Detect sources with npixels connected pixels at/above threshold in data smoothed by kernel
        # https://photutils.readthedocs.io/en/stable/segmentation.html
        self.segm_detect = detect_sources(self.data, detection_threshold, npixels=npixels, kernel=self.smooth_kernel)

        # Deblend: separate connected/overlapping sources
        # https://photutils.readthedocs.io/en/stable/segmentation.html#source-deblending
        self.segm_deblend = deblend_sources(self.data, self.segm_detect, npixels=npixels, kernel=self.smooth_kernel, nlevels=deblend_levels, contrast=deblend_contrast)

        if self.verbose:
            print(len(self.segm_detect.labels))


    def measure_source_properties(self, local_background_width=24, properties=output_properties):
        if version.parse(photutils.__version__) >= version.parse("1.4.0"):
            self.catalog = SourceCatalog(self.data-self.background_map.background, self.segm_deblend,
                                         convolved_data=self.convolved_data,  # photutils 1.4
                                         error=self.data_error, mask=self.data_mask,
                                         background=self.background_map.background, wcs=self.imwcs,
                                         localbkg_width=local_background_width)
        else:  # use filter_kernel instead of convolved_data
            self.catalog = SourceCatalog(self.data-self.background_map.background, self.segm_deblend,
                                         kernel=self.smooth_kernel,  # photutils < 1.4
                                         error=self.data_error, mask=self.data_mask,
                                         background=self.background_map.background, wcs=self.imwcs,
                                         localbkg_width=local_background_width)


        self.catalog_table = self.catalog.to_table(columns=properties)  # properties: quantities to keep

        # Convert fluxes to nJy units and to AB magnitudes
        for aperture in ['segment', 'kron']:
            flux    = self.catalog_table[aperture+'_flux']    * self.flux_units.to(u.nJy)
            fluxerr = self.catalog_table[aperture+'_fluxerr'] * self.flux_units.to(u.nJy)
            mag, magerr = fluxes2mags(flux, fluxerr)

            self.catalog_table[aperture+'_flux']    = flux
            self.catalog_table[aperture+'_fluxerr'] = fluxerr
            self.catalog_table[aperture+'_mag']     = mag
            self.catalog_table[aperture+'_magerr']  = magerr



    # def measure_background_and_rms(self):
    #
    #     self.measure_background_map()
    #
    #     self.background_rms = 1 / np.sqrt(self.wht)
    #
    #     if self.exptime:
    #
    #         self.exposure_time_map = self.exptime * self.background_map.background_rms_median**2 * self.wht
    #
    #         # effective gain parameter required to be positive everywhere (not zero), so adding small value 1e-8
    #         self.data_rms = calc_total_error(self.data, self.background_rms, self.exposure_time_map+1e-8)
    #
    #     else:
    #
    #         print('WARNING: Using naive data_rms calculation')
    #         self.data_rms = self.background_rms


    def get_rms(self):

        return self.data/self.data_rms

    def measure_rms(self):

        self.rms = self.get_rms()


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

    def make_cutout(self, x, y, width, extensions = ['sci', 'err', 'wht', 'bkg', 'bkg_rms']):
        """extract cut out"""

        if 'err' in extensions: err = np.zeros((width, width))
        if 'sci' in extensions: data = np.zeros((width, width))
        if 'wht' in extensions: wht = np.zeros((width, width))
        if 'bkg' in extensions: bkg = np.zeros((width, width))
        if 'bkg_rms' in extensions: bkg_rms = np.zeros((width, width))

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
        if 'err' in extensions: err[xstart:xend,ystart:yend] = self.err[xmin:xmax,ymin:ymax]
        if 'wht' in extensions: wht[xstart:xend,ystart:yend] = self.wht[xmin:xmax,ymin:ymax]
        if 'bkg' in extensions: bkg[xstart:xend,ystart:yend] = self.bkg[xmin:xmax,ymin:ymax]
        if 'bkg_rms' in extensions: bkg_rms[xstart:xend,ystart:yend] = self.bkg_rms[xmin:xmax,ymin:ymax]

        return ImageFromArrays(data, err = err, wht = wht, bkg = bkg, bkg_rms = bkg_rms, pixel_scale = self.pixel_scale, verbose = self.verbose)



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



    def flux(self, zeropoint = None):

        """ make flux image (not for photometry but for images) """

        if zeropoint:
            zp = zeropoint
        else:
            zp = self.zeropoint

        if zp:

            if self.units == 'e/s':

                return self.data/( 1E-9 * 10**(0.4*(zp-8.9)))

        else:
            print('WARNING: no zeropoint set')
            return


    def flux_rms(self, zeropoint = None):

        """ make flux image (not for photometry but for images) """

        if zeropoint:
            zp = zeropoint
        else:
            zp = self.zeropoint

        if zp:

            if self.units == 'e/s':

                return self.data_rms/( 1E-9 * 10**(0.4*(zp-8.9)))

        else:
            print('WARNING: no zeropoint set')
            return




class ImageFromMultiFITS(Image):

    def __init__(self, filename, idata = {'sci': 1, 'err': 2, 'wht': 4}, mask = None, verbose = False, mask_edge_thickness=10):

        """generate instance of image class from file"""

        self.verbose = verbose

        self.filename = filename
        self.hdu = fits.open(filename)
        self.header = self.hdu[0].header
        self.imwcs = wcs.WCS(self.hdu[idata['sci']].header, self.hdu)

        self.data = self.hdu[idata['sci']].data
        self.err = self.hdu[idata['err']].data
        self.wht = self.hdu[idata['wht']].data

        self.bkg = np.empty(self.data.shape)
        self.bkg_rms = np.empty(self.data.shape)

        # total error array (i.e., the background-only error plus Poisson noise due to individual sources)
        # https://photutils.readthedocs.io/en/stable/segmentation.html#photometric-errors
        #self.data_error = self.hdu[2].data  # 'ERR' i2d extension 2
        self.err = fits.open(self.filename)[idata['err']].data
        self.mask = np.isnan(self.err) # | (model.wht == 0)

        # Remove edge detections: Grow the mask by 10 pixels
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.binary_dilation.html
        self.mask = ndimage.binary_dilation(self.mask, iterations=mask_edge_thickness)

        # image_pixel_scale = np.abs(hdu[0].header['CD1_1']) * 3600
        self.pixel_scale = wcs.utils.proj_plane_pixel_scales(self.imwcs)[0]
        self.pixel_scale *= self.imwcs.wcs.cunit[0].to('arcsec')

        # --- need to generalise
        self.flux_units = JWST_flux_units * (self.pixel_scale * u.arcsec)**2

        if self.verbose:
            print(self.filename)
            ny, nx = self.data.shape
            outline = '%d x %d pixels' % (ny, nx)
            outline += ' = %g" x %g"' % (ny * self.pixel_scale, nx * self.pixel_scale)
            outline += ' (%.2f" / pixel)' % self.pixel_scale
            print(outline)





class ImageFromArrays(Image):

    def __init__(self, data, err = None, wht = None, bkg = None, bkg_rms = None, pixel_scale = None, verbose = True):

        """generate instance of image class from cutout"""

        print('here')

        self.verbose = verbose

        self.pixel_scale = pixel_scale
        self.data = data
        self.err = err
        self.wht = wht
        self.bkg = bkg
        self.bkg_rms = bkg_rms

        # for im in [self.data, self.wht, self.bkg, self.bkg_rms]:
        #     print(im.shape, np.median(im), np.mean(im), np.std(im), im.dtype)






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
