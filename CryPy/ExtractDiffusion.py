import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np
from CryPy.SignalProcessing import *
import scipy.ndimage
from lmfit.models import LorentzianModel, LinearModel
from lmfit import Parameters, Model

name = 'M78296NC5_Fluorescence_Trace'

from functools import partial


def lorentzian(f, f0, BG, FWHM, A):
    a = FWHM ** 2 / 2
    L = BG + A * a / (a + (f - f0) ** 2)
    return L


def fit_Line(f, y, model, params):
    model = Model(lorentzian)
    result = model.fit(y, params, f=f)
    return result.params


def fitted_Dif_map(data, mol='mol0', df=0.05, dt=1, molname='__', mean_BG='3000'):
    data['f'] = (data.Frequency - data.Frequency.mean()) / 1e9

    fmin = data.f.min()
    fmax = data.f.max()

    tmin = data.t.min()
    tmax = data.t.max()
    f = np.arange(fmin, fmax + df, df)
    # t = np.arange(tmin, tmax + dt, dt)
    xpix = int((fmax - fmin) // df) + 2
    ypix = int((tmax - tmin) // dt) + 2
    dmap = np.zeros([xpix, ypix], dtype=float)
    print(xpix, 'x', ypix)

    model = Model(lorentzian)
    params = Parameters()
    params.add('A', value=1000)
    params.add("f0", value=3)
    params.add("FWHM", value=2, min=0.5, max=3)
    params.add('BG', value=mean_BG, min=0, max=20000)

    for n in range(ypix - 1):
        T = (data.t > dt * n) & (data.t < dt * n + dt)
        ldata = data[T]
        newParams = fit_Line(ldata.f, ldata.mol0, model, params)
        params = newParams.copy()
        newParams['A'].value = 1  # normalizing each line. the intensity does not matter this way
        newParams['BG'].value = 0
        dmap[:, n] = model.eval(newParams, f=f)

    f, t = np.mgrid[slice(fmin, fmax + df, df), slice(tmin, tmax + dt, dt)]

    fig, axs = plt.subplots(1, 1)

    ax = axs
    c = ax.pcolor(f, t, dmap, cmap='viridis', vmin=0, vmax=1)
    ax.set_title('Spectral Diffusion')
    ax.set_xlabel('$\Delta f$  (GHz)')
    ax.set_ylabel('Time (s)')
    fig.colorbar(c, ax=ax)
    plt.show()


def Extract_diffusion(data, mol='mol0', df=0.5, dt=1, molname='__', mean_BG='3000'):
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

    print("Image is {}x{}".format(xpix, ypix))
    z = np.zeros([xpix, ypix], dtype=float) + mean_BG
    for ind in data.index:
        tpix = int((data.ElapsedTime[ind]) // dt) - 1
        fpix = int((data.Frequency[ind] - fmin) // df) - 1
        z[fpix, tpix] = np.max([data[mol][ind], z[fpix, tpix]])

    z_min, z_max = np.min(z), np.abs(z).max()
    # plt.imshow(z)
    fig, axs = plt.subplots(1, 2)

    # z = z[:-2, :-2]
    ax = axs[0]
    # z = scipy.ndimage.zoom(z,3)
    # f = scipy.ndimage.zoom(f, 3)
    # t = scipy.ndimage.zoom(t, 3)
    c = ax.pcolor(f, t, z, cmap='inferno', vmin=2000, vmax=z_max)
    ax.set_title('Spectral Diffusion')
    ax.set_xlabel('$\Delta f$  (GHz)')
    ax.set_ylabel('Time (s)')
    fig.colorbar(c, ax=ax)

    zf = butter_lowpass_filter(z[:], 0.7, order=5)
    ax = axs[1]
    # z= z[:-2, :-2]
    # zf = zf[:-2, :-2]
    c = ax.contourf(f, t, zf, cmap='inferno', vmin=z_min, vmax=z_max)
    ax.set_title('LowPass Filtered in each line')
    ax.set_title('Spectral Diffusion')
    ax.set_xlabel('$\Delta f$  (GHz)')
    ax.set_ylabel('Time (s)')
    fig.colorbar(c, ax=ax)

    fig.tight_layout()
    fig.savefig(molname + '_DataOutPut\\' + molname + r'_SpectralDiffusion.png')
    plt.show()
