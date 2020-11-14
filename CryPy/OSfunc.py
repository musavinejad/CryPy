import os
import glob
import pandas as pd

def prepare_folders(MolNameTag):
    print('Preparing the folders for the results.')
    folder = os.path.abspath('main.py')[:-7]
    os.chdir(folder)
    print('Current System Path Changed to : ' + folder)
    os.makedirs(MolNameTag + "_DataOutPut", exist_ok=True)
    print('Folder <' +MolNameTag + "_DataOutPut> created.")
    return folder


def select_data():
    files = pd.DataFrame()
    ff = glob.glob("Data/*.cam")
    files['name'] = ff
    print(files)
    print('Select your data by index :')
    a = int(input())
    name = files.name[a]
    name = name[:-11]  ## taking out "_imagem.cam"
    print('analysing ', name[5:])
    return name

