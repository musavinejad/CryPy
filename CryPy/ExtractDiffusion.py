import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np


name = 'Molecules_Fluorescence_Trace'


def Extract_diffusion(data, mol = 'mol1',clockRate= 50):
    print('Hi im diffusion mapper!')
    #plt.plot(data.Frequency/1e9,data[mol])
    fig, ax = plt.subplots(1)
    f = np.array([data.Frequency[-1000:-1] ,data.Frequency[-1001:]])
    t = np.array([data.ElapsedTime[-1000:-1], data.ElapsedTime[-1001:]])
    f,t = np.mgrid[f,t]
    fluo = np.zeros(f.shape)
    plt.pcolor(f,t, data[mol] ,cmap='RdBu')











print('Hi!')
folder = os.path.abspath(__file__)[:-(('ExtractDiffusion.py').__len__())]
os.chdir(folder)
os.makedirs("Spectral_Diffusion_Map", exist_ok=True)
data= pd.read_csv(name + '.csv', delimiter=',')

print('analysing ', name[5:])

Extract_diffusion(data, 'mol8')
