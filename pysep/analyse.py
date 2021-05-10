

import os
import sys
import h5py
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import FLARE.plt as fplt



def default_plot():

    fig = plt.figure(figsize = (4,4))

    left  = 0.15
    bottom = 0.15
    height = 0.8
    width = 0.8

    ax = fig.add_axes((left, bottom, width, height))

    return fig, ax



class Analyser:

    def __init__(self, cat_file, show_plots = True, save_plots = False):

        self.show_plots = show_plots
        self.save_plots = save_plots

        self.cat = h5py.File(cat_file, 'r')

        # self.plot_dir = f'{output_dir}/{output_filename}'
        #
        # if not os.path.exists(self.plot_dir):
        #     os.makedirs(self.plot_dir)


    def explore_hdf5(self):

        def get_name_shape(name, item):
            shape = ''
            if hasattr(item, 'value'):
                shape = item.shape
            print(name, shape)

        self.cat.visititems(get_name_shape)


    def scatter(self, x_name, y_name, z_name = None, logx = False, logy = False, logz = False, s = 10, cmap = 'inferno'):

        fig, ax = default_plot()

        if logx:
            x = np.log10(self.cat[x_name])
        else:
            x = self.cat[x_name]

        if logy:
            y = np.log10(self.cat[y_name])
        else:
            y = self.cat[y_name]

        if z_name:

            if logz:
                z = np.log10(self.cat[z_name])
            else:
                z = self.cat[z_name]

            norm = mpl.colors.Normalize(vmin = np.min(z), vmax = np.max(z))

            ax.scatter(x, y, c = norm(z), cmap = 'inferno',  s = s)

        else:

            ax.scatter(x, y, s = s)

        ax.set_xlabel(x_name)
        ax.set_ylabel(y_name)

        if self.show_plots:
            plt.show()
