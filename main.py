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

molname = 'M76876NC5'

#########################################

############     Local codes constants:

name = ""
neighborhood_size = 5
threshold = 1500

############


# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    folder= prepare_folders(MolNameTag = molname)

    name = select_data()  ## put the data you want to analyze in a folder named "...\Data\"
    reader = read_camera(location=name + '_imagem.cam', load_all=False)
    #reader = np.array(read_camera(location=name + '_imagem.cam', load_all=True))  # read data as np.array

    data= pd.read_csv(name + '.dat', delimiter='\t', usecols=['DAQ_Index','ElapsedTime','Frequency'])
    mean_BG, STD_BG = measure_noise(reader)  # for the threshold in finding molecules
    max_img = get_Max_Trace(reader)
    ### plotting the max trace image and the molecules found in it.
    fig, ax = plt.subplots(1)
    ax.imshow(max_img, vmin= mean_BG, vmax= mean_BG+10* STD_BG )
    fig.suptitle(molname+ '_Max Trace')

    fig.savefig(molname + '_DataOutPut\\' + molname + r'MaxTrace.png' )
    molx, moly = find_mols(max_img, neighborhood_size=5, threshold=4* STD_BG)
    #molx[1] = 46
    #moly[1] = 81
    #molx = np.array([82])
    #moly = np.array([70])
    add_patches_to_img(ax, molx, moly)  ## Adds rectangles around the molecules
    fig.savefig(molname + '_DataOutPut\\' + molname +r'Maxtrace_withMolecules.png')

    ### max trace image is saved. ,molx/y contains the coordinates of possible molecules

    reader = np.array(reader.read_all())

    get_mols_spectrum(data, reader, molx, moly , molname)

    #data=data.sort_values(by=['Frequency'])
    plot_mols(data, molx.size, molname)

    fig, ax = plt.subplots(1)
    plt.xlabel('Detuning (GHz)')
    plt.ylabel('Fluorescence (A.U.)')
    ax.plot((data.Frequency-data.Frequency.mean())/1e9 , data.mol8/mean_BG + (data.ElapsedTime//0.5))
    #Extract_diffusion(data, 'mol0' , df= 1.5,dt= 0.3, molname= molname )

    #butter_lowpass_filter(z[1], cutoff, fs, order = 5)






