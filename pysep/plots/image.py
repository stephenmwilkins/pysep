


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm



def make_flux_plot(img, vmin = None, vmax = None, ext = 'sci', scaling = False, cmap = cm.magma):


    fig, ax = plt.subplots(1, 1, figsize = (4,4))
    plt.subplots_adjust(left=0, top=1, bottom=0, right=1, wspace=0.01, hspace=0.0)

    if ext == 'sci': im_ = img.sci
    if ext == 'wht': im_ = img.wht

    if scaling:
        im = scaling(im_)
    else:
        im = im_ # linear scaling

    if not vmin:
        vmin = np.min(im)

    if not vmax:
        vmax = np.max(im)

    ax.set_axis_off()
    ax.imshow(im, cmap = cmap, vmin = vmin, vmax = vmax, origin = 'lower') # choose better scaling

    plt.show()
    plt.close(fig)



def make_flux_plots(imgs, vmin = 0, vmax = None, show = True):

    n = len(imgs)

    fig, axes = plt.subplots(1, n, figsize = (4*n,4))
    plt.subplots_adjust(left=0, top=1, bottom=0, right=1, wspace=0.01, hspace=0.0)

    filters = list(imgs.keys())

    if not vmax:
        vmax = np.max([np.max(imgs[f].sci) for f in filters])

    for ax, f in zip(axes, filters):
        img = imgs[f]
        ax.set_axis_off()
        ax.imshow(img.sci, cmap = cm.magma, vmin = vmin, vmax = vmax, origin = 'lower')

    plt.show()
    plt.close(fig)


# ---------------------------------------------------------------------------------------------

def make_significance_panel(ax, img, threshold = 2.5):

    sig = (img.sci/img.noise)

    ax.imshow(sig, cmap = cm.Greys, vmin = -5.0, vmax = 5.0, origin = 'lower')
    ax.imshow(np.ma.masked_where(sig <= threshold, sig), cmap = cm.plasma, vmin = threshold, vmax = 100, origin = 'lower')
    ax.set_axis_off()

    return ax


def make_significance_plot(img, patches = None, threshold = 2.5, save_file = None, show = True):

    fig, ax = plt.subplots(1, 1, figsize = (4,4))
    plt.subplots_adjust(left=0, top=1, bottom=0, right=1, wspace=0.01, hspace=0.0)

    ax = make_significance_panel(ax, img, threshold = threshold)

    if patches:
        for patch in patches:
            ax.add_artist(patch)
            patch.set_alpha(0.5)

    plt.show()
    plt.close(fig)


def make_significance_plots(imgs, threshold = 2.5, save_file = None, show = True):

    n = len(imgs)

    fig, axes = plt.subplots(1, n, figsize = (4*n,4))
    plt.subplots_adjust(left=0, top=1, bottom=0, right=1, wspace=0.01, hspace=0.0)

    filters = list(imgs.keys())

    for ax, f in zip(axes, filters):
        ax = make_significance_panel(ax, imgs[f])

    if save_file: fig.savefig(f'{save_file}')
    if show: plt.show()





# ---------------------------------------------------------------------------------------------

def make_segm_panel(ax, segm):

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

    ax = make_segm_panel(ax, segm)

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
