import numpy as np
from scipy.signal import butter,filtfilt
from scipy.signal import savgol_filter

def butter_lowpass_filter(z, cutoff, order = 5):
    # Get the filter coefficients
    print(z)
    b, a = butter(order, cutoff, btype='low', analog=False)
    y = filtfilt(b, a, z, axis=0)
    return y


def add_smothed_freq(data):
    # data['rawFreq'] = data['Frequency']
    data['f'] = (data['Frequency']-data['Frequency'].mean())/1e9
    data['f'] = savgol_filter(data['f'], 3, 1)

def add_smothed_time(data):
    # data['rawFreq'] = data['Frequency']
    data['t'] = savgol_filter(data['ElapsedTime'], 3, 1)