import numpy as np
import matplotlib.pyplot as mpl
import matplotlib.animation as ani
from pylablib.aux_libs.file_formats import waveguide, cam
from pylablib.core.datatable import table
from pylablib.core.utils import files as file_utils, dictionary, plotting
from pylablib.core.fileio import loadfile, savefile
import os.path
import glob
import pandas as pd
import scipy
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
import matplotlib.pyplot as plt
import matplotlib.patches as patches

#############################
def select_data():
    files = pd.DataFrame()
    ff = glob.glob("Data/*.cam")
    files['name'] = ff
    print(files)
    print('Select your data by index :')
    a = int(input())
    return files.name[a]


#############################
def read_camera(location, load_all=False):
    reader = cam.CamReader(location)
    if load_all:
        reader = reader.read_all()
    return reader


#############################
def get_Max_Trace(reader):
    return np.max(reader, axis=0)


#############################
def find_mols(data, neighborhood_size, threshold):
    data_max = filters.maximum_filter(data, neighborhood_size)
    maxima = (data == data_max)
    data_min = filters.minimum_filter(data, neighborhood_size)
    diff = ((data_max - data_min) > threshold)
    maxima[diff == 0] = 0

    labeled, num_objects = ndimage.label(maxima)
    slices = ndimage.find_objects(labeled)
    x, y = [], []
    for dy, dx in slices:
        x_center = (dx.start + dx.stop - 1) / 2
        x.append(x_center)
        y_center = (dy.start + dy.stop - 1) / 2
        y.append(y_center)
    return np.array(x, dtype=np.uint8), np.array(y,dtype=np.uint8)



#############################
def measure_noise(reader, x=20, y=20):
    mean = np.mean(reader[100:1000, x, y])
    print("mean value = %.2f" %(mean) )
    std = np.std(reader[100:1000, x, y])
    print("Background STD = %.2f" % (std))
    return float(mean), float(std)

#############################
def save_figure(img, figttl):
    fig,ax = plt.subplot(1)

#############################
def add_patches_to_img(ax, molx, moly):
    for i in range(0, molx.size):
        ax.add_patch(
            patches.Rectangle([molx[i] - 2, moly[i] - 2], 4, 4, fill=False, edgecolor='white', linestyle='--'))
        plt.text(molx[i] - 2, moly[i] - 2, '{}'.format(i), color='orange')


#############################
def plot_mols_spectrum(data, reader, molx, moly):
    n = molx.size
    print('Plotting {} molecules spectrum...'.format(n))
    for i in range (0,n):
        fluor_trace = np.mean(np.mean(reader[:, moly[i]-2 : moly[i]+2, molx[i]-2 : molx[i]+2], axis=1), axis=1)
        data['mol'+str(i)] = fluor_trace


def plot_tera(new_data,reader,center=(5,5),rngs=(3,3),save_name="tera_scan"):
    folder = "Tera_scans"
    # Create target Directory if don't exist
    if not os.path.exists(folder):
        os.mkdir(folder)
        print("Directory ", folder, " Created ")
    # create the path, where the video gets saved
    path = os.path.join(folder, save_name+".png")

    data=new_data.copy()
    frames=np.array(reader.read_all())
    fluor_trace=np.mean(np.mean(frames[:,center[0]-rngs[0]//2:center[0]+rngs[0]//2,center[1]-rngs[1]//2:center[1]+rngs[1]//2],axis=1),axis=1)
    data.c.append(fluor_trace,names="Fluor") #appends a new colum to data, new colum data is fluor_trace, name is "Fluor"
    waveguide.trim_jumps(data,jump_size=500*1e6,trim=3,x_column="Frequency")
    data=data.sort_by("Frequency")

    f=mpl.figure(figsize=(14,6))
    ax=f.add_subplot(111)
    ax.plot(data["Frequency"]/1e12,data["Fluor"])
    ax.set_xlim([382, 382.38])
    ax.set_xlabel("Frequency [THz]")
    ax.set_ylabel("Fluorescence [A. U.]")
    f.savefig(path,dpi=300)
    np.savetxt(save_name+"raw.txt" , data , delimiter=',')
    mpl.close(f)