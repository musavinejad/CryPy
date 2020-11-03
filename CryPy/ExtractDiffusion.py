import matplotlib.pyplot as plt
import pandas as pd
import os

name = 'Molecules_Fluorescence_Trace'


def Extract_diffusion(data, mol = 'mol1',clockRate= 50):
    print('Hi im diffusion mapper!')
    #plt.plot(data.Frequency/1e9,data[mol])
    data.ElapsedTime.max











print('Hi!')
folder = os.path.abspath(__file__)[:-(('ExtractDiffusion.py').__len__())]
os.chdir(folder)
os.makedirs("Spectral_Diffusion_Map", exist_ok=True)
data= pd.read_csv(name + '.csv', delimiter=',')

print('analysing ', name[5:])

Extract_diffusion(data, 'mol8')
