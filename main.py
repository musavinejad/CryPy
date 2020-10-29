import os.path
from typing import Optional, List, Any, Union

import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from pylablib.aux_libs.file_formats.cam import CamReader
from scipy.signal import find_peaks
from pylablib.aux_libs.file_formats import waveguide, cam
from pylablib.core.fileio import loadfile, savefile
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
from CryPy import *

import matplotlib.patches as patches

############     constants:

name = ""
neighborhood_size = 5
threshold = 1500

############


# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    folder = os.path.abspath(__file__)[:-7]
    os.chdir(folder)
    os.makedirs("DataOutPut", exist_ok=True)

    name = select_data()
    name = name[:-11]  # taking out "_imagem.cam"
    print('analysing ', name[5:])
    reader = np.array(read_camera(location=name + '_imagem.cam', load_all=True))  # read data as np.array
    mean_BG, STD_BG = measure_noise(reader, 40, 40)  # for the threshold
    max_img = get_Max_Trace(reader)
        ### plotting the max trace image and the molecules found in it.
    fig, ax = plt.subplots(1)
    ax.imshow(max_img)
    os.chdir(folder + r'\DataOutput')
    fig.savefig(r'Maxtrace.png')
    molx, moly = find_mols(max_img, neighborhood_size=5, threshold=mean_BG + 3 * STD_BG)
    for i in range(0, molx.size):
        ax.add_patch(
            patches.Rectangle([molx[0, i] - 2, moly[0, i] - 2], 4, 4, fill=False, edgecolor='white', linestyle='--'))
        plt.text(molx[0, i] - 2, moly[0, i] - 2, '{}'.format(i), color='orange')
    fig.savefig(r'Maxtrace_withMolecules.png')
        ### max trace image is saved. ,molx/y contains the coordinates of possible molecules