import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np
from CryPy.SignalProcessing import *

name = 'M78296NC5_Fluorescence_Trace'


def Extract_diffusion(data, mol = 'mol0', df=0.5, dt =1 , molname= '__'):
    print('mapping the spectral diffusion')

    # generate 2 2d grids for the x & y bounds
    data.Frequency = (data.Frequency - data.Frequency.mean()) / 1e9
    fmin = data.Frequency.min()
    fmax = data.Frequency.max()

    tmin = data.ElapsedTime.min()
    tmax = data.ElapsedTime.max()

    f, t = np.mgrid[slice(fmin, fmax + df, df), slice(tmin, tmax + dt, dt)]
    xpix = int((fmax - fmin) // df) + 2
    ypix = int((tmax - tmin) // dt) + 2

    z = np.zeros([xpix, ypix], dtype=float)
    for ind in data.index:
        tpix = int((data.ElapsedTime[ind]) // dt) - 1
        fpix = int((data.Frequency[ind] - fmin) // df) - 1
        z[fpix, tpix] = data[mol][ind]


    z_min, z_max = np.min(z), np.abs(z).max()

    fig, axs = plt.subplots(1)



    z = z[:-1, :-1]
    ax = axs
    c = ax.pcolor(f, t, z, cmap='inferno', vmin=2000, vmax=z_max)
    ax.set_title('Spectral Diffusion')
    ax.set_xlabel('$\Delta f$  (GHz)')
    ax.set_ylabel('Time (s)')
    fig.colorbar(c, ax=ax)


    # zf= butter_lowpass_filter(z[:], 0.9, order = 5)
    # ax = axs[1]
    # # z= z[:-2, :-2]
    # c = ax.pcolormesh(f, t, zf, cmap='inferno', vmin=z_min, vmax=z_max)
    # ax.set_title('LowPass Filtered in each line')
    # fig.colorbar(c, ax=ax)

    fig.tight_layout()
    fig.savefig(molname + '_DataOutPut\\' + molname + r'_SpectralDiffusion.png')
    plt.show()
    #
    # fig, axs = plt.subplots(2)
    # ax = axs[0]
    # ax.plot(t[1,:-1],z[1])
    # ax.plot(t[1,:-1],z[10])
    #
    # # ax = axs[1]
    # # ax.plot(f[:-1,1], z[:,1])
    # # ax.plot(f[:-1,1], z[:,10])


