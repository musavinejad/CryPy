import os.path
from typing import Optional, List, Any, Union

import pandas as pd
from docutils.nodes import figure
from matplotlib import pyplot as plt
import numpy as np
from pylablib.aux_libs.file_formats.cam import CamReader
from scipy.signal import find_peaks
from pylablib.aux_libs.file_formats import waveguide, cam
from pylablib.core.fileio import loadfile, savefile
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
from CryPy import *
import matplotlib
from CryPy.OSfunc import *
from CryPy.ExtractDiffusion import *
from CryPy.SignalProcessing import *

########################################
## input documentation names :

molname = 'M76824NC5_1V_3'

#########################################

############     Local codes constants:

name = ""
neighborhood_size = 5
threshold = 1500


############


class ScData:

    # name = ''
    # reader = read_camera
    # data = pd.DataFrame
    # mean_BG = 0
    # STD_BG = 0
    # max_img = np.array
    # molx, moly = np.array([]) , np.array([])
    # folder = ''
    def __init__(self, molname="", Auto = False):
        self.molname = molname
        self.read()
        if Auto:
            self.find_molecules()
            self.add_mol_spec()

    def read(self):
        self.folder = prepare_folders(MolNameTag=self.molname)
        self.name = select_data()  ## put the data you want to analyze in a folder named "...\Data\"
        self.reader = read_camera(location=self.name + '_imagem.cam',
                                  load_all=False)

        self.data = pd.read_csv(self.name + '.dat', delimiter='\t',
                                usecols=['DAQ_Index', 'ElapsedTime', 'Frequency'])
        add_smothed_freq(self.data)
        add_smothed_time(self.data)
        self.mean_BG, self.STD_BG = measure_noise(self.reader)
        self.max_img = get_Max_Trace(self.reader)

    def find_molecules(self , threshold= 5000):
        self.molx, self.moly = find_mols(self.max_img, neighborhood_size=5,
                                         threshold=threshold)

    def add_mol_spec(self):
        get_mols_spectrum(self.data, np.array(self.reader.read_all()),
                          self.molx, self.moly, self.molname)

    def plot_Diffusion(self, molecule='mol0', df=0.5, dt=10):
        Extract_diffusion(self.data.copy(), molecule, df, dt,
                          molname=self.molname, mean_BG=self.mean_BG)

    def plot_mols(self):
        plot_mols(self.data, self.molx.size, self.molname)

    def plot(self):
        self.mximg_fig, self.mximg_ax = plt.subplots(1)
        self.mximg_ax.imshow(self.max_img, vmin=self.mean_BG, vmax=3 * self.mean_BG)
        self.mximg_fig.suptitle(self.molname + '_Max Trace')

        self.mximg_fig.savefig(self.molname + '_DataOutPut\\' + self.molname + r'MaxTrace.png')
        add_patches_to_img(self.mximg_ax, self.molx, self.moly)  ## Adds rectangles around the molecules
        # fig.savefig(molname + '_DataOutPut\\' + molname +r'Maxtrace_withMolecules.png')

    def take_molecules(self):
        self.fig, self.ax = plt.subplots(1)
        self.ax.imshow(self.max_img, vmin=self.mean_BG, vmax=3 * self.mean_BG)
        self.fig.suptitle(self.molname + '_Max Trace')
        plt.show()
        plt.waitforbuttonpress()

        # while True:
        pts = []
        plt.title('Select the molecules with the mouse:')
        pts = (plt.ginput(n=-1, show_clicks=True, mouse_stop=3, timeout=-1))

        # ph = plt.fill(pts[:, 0], pts[:, 1], 'r', lw=2)
        pts2 = [list(i) for i in pts]
        # print('Press Enter to stop adding.\n', self.pts2)

        self.molx = np.array([round(x) for x in np.asarray(pts2)[:, 0]], dtype= int)
        self.moly = np.array([round(x) for x in np.asarray(pts2)[:, 1]], dtype= int)
        print(self.molx, self.moly)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    # matplotlib.use("Qt5Agg")
    a = ScData(molname="ss" , Auto= True)
    data = a.data.copy()

    # a.take_molecules()

