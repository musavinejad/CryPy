
import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np
from CryPy.SignalProcessing import *
import scipy.ndimage
from lmfit.models import LorentzianModel, LinearModel
from lmfit import Parameters, Model

name = 'M78296NC5_Fluorescence_Trace'


cdef float lorentzian(f, f0, BG, FWHM, A):
    a = FWHM ** 2 / 2
    L = BG + A * a / (a + (f - f0) ** 2)
    return L


cdef fit_Line(f, y, model, params):
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
    xpix = int((fmax - fmin) // df) +2
    ypix = int((tmax - tmin) // dt) +2
    dmap = np.zeros([xpix, ypix], dtype=float)
    print(xpix, 'x', ypix)

    model = Model(lorentzian)
    params = Parameters()
    params.add('A', value=1000)
    params.add("f0", value=3)
    params.add("FWHM", value=2, min=0.5, max=3)
    params.add('BG', value=int(mean_BG), min=0, max=20000)

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
