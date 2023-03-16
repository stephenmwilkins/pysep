


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm



def make_mutlipanel_image(n, size = 3):

    fig, axes = plt.subplots(1, n, figsize = (size*n, size))
    plt.subplots_adjust(left=0, top=1, bottom=0, right=1, wspace=0.01, hspace=0.0)

    return fig, axes



def img_panel(ax, im, vmin = None, vmax = None, scaling = False, cmap = cm.magma):

    if not vmin:
        vmin = np.min(im)
    if not vmax:
        vmax = np.max(im)

    if scaling:
        im = scaling(im)

    ax.set_axis_off()
    ax.imshow(im, cmap = cmap, vmin = vmin, vmax = vmax, origin = 'lower') # choose better scaling

    return ax


def make_image_plot(im, vmin = None, vmax = None, scaling = False, cmap = cm.magma):

    fig, ax = plt.subplots(1, 1, figsize = (4,4))
    plt.subplots_adjust(left=0, top=1, bottom=0, right=1, wspace=0.01, hspace=0.0)

    ax = img_panel(ax, im, vmin = vmin, vmax = vmax, scaling = scaling, cmap = cmap)

    return fig, ax

def make_images_plot(ims, vmin = None, vmax = None, scaling = False, cmap = cm.magma):

    n = len(ims)
    fig, axes = plt.subplots(1, n, figsize = (4*n,4))
    plt.subplots_adjust(left=0, top=1, bottom=0, right=1, wspace=0.01, hspace=0.0)

    for im, ax in zip(ims, axes):
        ax = img_panel(ax, im, vmin = vmin, vmax = vmax, scaling = scaling, cmap = cmap)

    return fig, axes




def make_flux_plots(imgs, vmin = 0, vmax = None, show = True):

    n = len(imgs)

    fig, axes = plt.subplots(1, n, figsize = (4*n,4))
    plt.subplots_adjust(left=0, top=1, bottom=0, right=1, wspace=0.01, hspace=0.0)

    filters = list(imgs.keys())

    if not vmax:
        vmax = np.max([np.max(imgs[f].data) for f in filters])

    for ax, f in zip(axes, filters):
        img = imgs[f]
        ax.set_axis_off()
        ax.imshow(img.data, cmap = cm.magma, vmin = vmin, vmax = vmax, origin = 'lower')

    plt.show()
    plt.close(fig)


# ---------------------------------------------------------------------------------------------







def significance_panel(ax, img, threshold = 2.5):

    sig = (img.data-img.bkg)/img.bkg_rms
    # sig = (img.data)/img.err

    ax.imshow(sig, cmap = cm.Greys, vmin = -threshold*2, vmax = threshold*2, origin = 'lower')
    ax.imshow(np.ma.masked_where(sig <= threshold, sig), cmap = cm.plasma, vmin = threshold, vmax = 100, origin = 'lower')
    ax.set_axis_off()

    return ax


def make_significance_plot(img, patches = None, threshold = 2.5, save_file = None, show = True):

    fig, ax = plt.subplots(1, 1, figsize = (4,4))
    plt.subplots_adjust(left=0, top=1, bottom=0, right=1, wspace=0.01, hspace=0.0)

    ax = significance_panel(ax, img, threshold = threshold)

    if patches:
        for patch in patches:
            ax.add_artist(patch)
            patch.set_alpha(0.5)

    return fig, ax


def make_significance_plots(imgs, threshold = 2.5, save_file = None, show = True):

    n = len(imgs)

    fig, axes = plt.subplots(1, n, figsize = (4*n,4))
    plt.subplots_adjust(left=0, top=1, bottom=0, right=1, wspace=0.01, hspace=0.0)

    filters = list(imgs.keys())

    for ax, f in zip(axes, filters):
        ax = significance_panel(ax, imgs[f])

    return fig, ax





# ---------------------------------------------------------------------------------------------

def segm_panel(ax, segm):

    # --- create random colour map
    vals = np.linspace(0, 1, 256)
    np.random.shuffle(vals)
    random_cmap = plt.cm.colors.ListedColormap(plt.cm.jet(vals))

    ax.imshow(segm.data, cmap = random_cmap, origin = 'lower')
    ax.set_axis_off()

    return ax


def make_segm_plot(segm, imsize = 4, save_file = None, show = True):

    fig, ax = plt.subplots(1, 1, figsize = (imsize,imsize))

    plt.subplots_adjust(left=0, top=1, bottom=0, right=1, wspace=0.0, hspace=0.0)

    ax = segm_panel(ax, segm)

    if show:
        plt.show()

    plt.close(fig)



def make_significance_segm_plot(img, segm, threshold = 2.5, save_file = None, show = True):


    fig, axes = plt.subplots(1, 2, figsize = (8,4))
    plt.subplots_adjust(left=0, top=1, bottom=0, right=1, wspace=0.01, hspace=0.0)

    ax_sig, ax_segm = axes

    ax_sig = make_significance_panel(ax_sig, img)
    ax_segm = make_segm_panel(ax_segm, segm)


    if show:
        plt.show()
    plt.close(fig)




def make_rgb_image(r, g, b, Q=5, stretch=0.02, show = True, save = False, imsize = 3):


    from astropy.visualization import make_lupton_rgb

    fig, ax = plt.subplots(1, 1, figsize = (imsize,imsize))

    plt.subplots_adjust(left=0, top=1, bottom=0, right=1, wspace=0.0, hspace=0.0)

    # rgb_image = make_lupton_rgb(r, g, b, Q=5, stretch=0.02)
    # rgb_image = make_lupton_rgb(r, g, b)

    ax.imshow(np.array([r,g,b]).T)

    if show:
        plt.show()

    plt.close(fig)
