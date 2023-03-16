

import os
import sys

import numpy as np

import pysep.sep as sep
import pysep.utils
import pysep.plots.image
import pysep.analyse



# --- read the test data. The test data is organised as: data/test/{f}_{image_type}.fits where f is the filter and image_type is e.g. sci, wht

filters = ['f435w','f606w','f775w','f814w','f850lp','f105w','f125w','f140w','f160w']
filters = ['f435w','f814w','f125w']

zeropoints = {'f105w': 26.269, 'f125w': 26.230, 'f140w': 26.452, 'f160w': 25.946, 'f435w': 25.684, 'f606w': 26.505, 'f775w': 25.678, 'f814w': 25.959, 'f850lp': 24.867}


# data_dir = 'data/test'
# filename = lambda f, ext: f'{data_dir}/{f}_{ext}' # define a function that gives the full filename for different filters
# imgs = {f: pysep.utils.ImageFromFITS(filename(f,'sci'), filename(f,'wht'), calculate_background = False, zeropoint = zeropoints[f]) for f in filters}
#

data_dir = '/Users/stephenwilkins/Dropbox/Research/data/images/hubble/xdf'
sci_filename = lambda f: f'{data_dir}/{f}_sci_convolved'
wht_filename = lambda f: f'{data_dir}/{f}_wht'
imgs = {f: pysep.utils.ImageFromFITS(sci_filename(f), wht_filename(f), zeropoint = zeropoints[f]) for f in filters}


# r = pysep.utils.create_stack({f:imgs[f] for f in ['f105w','f125w','f140w','f160w']})
# g = pysep.utils.create_stack({f:imgs[f] for f in ['f775w','f814w','f850lp']})
# b = pysep.utils.create_stack({f:imgs[f] for f in ['f435w','f606w']})


r = imgs['f125w'].flux()
g = imgs['f814w'].flux()
b = imgs['f435w'].flux()



maxf = 10

for im, sf in zip([r,g,b],[0.5, 0.5, 1.]):
    print(np.max(im))
    im *= sf
    im[im>maxf] = maxf
    im /= maxf




# for im in [r,g,b]:
#     sd = np.std(im.data[np.fabs(im.sn())<5])
#     im.data[im.data>sd*10] = sd*10
#     print(np.max(im.data))


#
# r.data /= np.sum(r.data)
# g.data /= np.sum(g.data)
# b.data /= np.sum(b.data)
#
# r.data /= np.max(r.data)
# g.data /= np.max(g.data)
# b.data /= np.max(b.data)

# r.data[r.data!=r.data] = 0.
# g.data[g.data!=g.data] = 0.
# b.data[b.data!=b.data] = 0.
#

# rgb = r + g + b
#
# r /= np.sum(r)
# g /= np.sum(g)
# b /= np.sum(b)
#
# r /= np.sum(rgb)
# g /= np.sum(rgb)
# b /= np.sum(rgb)
#
# r /= np.max(r)
# g /= np.max(g)
# b /= np.max(b)
#
# r = np.log10(r)
# g = np.log10(g)
# b = np.log10(b)

from astropy.visualization import make_lupton_rgb
import matplotlib.pyplot as plt


imsize = 1
fig, ax = plt.subplots(1, 1, figsize = (imsize,imsize))
plt.subplots_adjust(left=0, top=1, bottom=0, right=1, wspace=0.0, hspace=0.0)

# rgb_image = make_lupton_rgb(r, g, b, Q=10, stretch=0.02)
# ax.imshow(rgb_image)

ax.imshow(np.array([r,g,b]).T)
ax.set_axis_off()
# ax.imshow(np.array([r.data,g.data,b.data]).T)

plt.savefig('rgb.png', dpi = r.shape[0])
