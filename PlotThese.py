import os
import pandas as pd
from matplotlib import pyplot as plt
filelist = os.listdir()

path = os.path.abspath(__file__)
os.makedirs(path[:-2] +"/Plots/" , exist_ok= True)
for x in filelist:
    if x.endswith("settings.dat"):
        continue
    if x.endswith(".dat"):
        data = pd.read_csv(x, delimiter="\t")
        plt.ylabel('Counts/S')
        plt.xlabel('THz')
        fig = plt.plot(data.Frequency, data.APD)
        plt.savefig('Plots/'+x[:-4]+'.png')
        plt.close()


for x in filelist:
    if x.endswith("settings.dat"):
        continue
    if x.endswith(".dat"):
        data = pd.read_csv(x, delimiter="\t")
        fig = plt.plot(data.Frequency, data.APD)
plt.savefig('Plots/'+'TotalScan.png')

