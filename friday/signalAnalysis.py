## IMPORTS
import numpy as np
import datetime
from os.path import isdir, getsize
from glob import glob
from sys import exit
import matplotlib.pylab as plt
import math
import DCIPtools as DCIP
from scipy import fftpack
import spectrum

def synthetic_current(sample_freq, time_base):
    # Creation of synthetic current signal of amplitude 1.
    # Input:
    # sample_freq: desired sample frequency of output current file
    # time_base: desired time base of synthetic current signal
    # Output:
    # time: time array
    # syn_current: time serie of synthetic signal
    N = 16384
    time = np.linspace(1, N, N) / sample_freq
    N_cycle = time_base * 4
    syn_current = np.asarray([0. for i in range(N)])
    for i in range(N):
        if np.mod(time[i], N_cycle) < time_base:
            syn_current[i] = 1.
        elif time_base <= np.mod(time[i], N_cycle) < time_base * 2:
            syn_current[i] = 0.
        elif time_base * 2 <= np.mod(time[i], N_cycle) < time_base * 3:
            syn_current[i] = -1.
        elif time_base * 3 <= np.mod(time[i], N_cycle) < time_base * 4:
            syn_current[i] = -0.

    return time, syn_current

def get_maxima(period, nfft, freq, signal_fft):
    # Get local maxima and minima of current waveform using the derivative
    # Creation of synthetic current signal of amplitude 1.
    # Input:
    # period: not used anymore, to be removed
    # nfft: not used anymore, to be removed
    # freq: frequency array corresponding to the input FFT
    # signal_fft: FFT of signal to get maxima
    # Output:
    # periods: frequency array corresponding to maxima locations

    derivative = [signal_fft[i + 1] - signal_fft[i] for i in range(int(len(freq) / 2))]

    periods = []
    index = []
    for i in range(int(len(freq) / 2) - 1):
        if (derivative[i] > 0 and derivative[i + 1]) < 0:
            periods.append(freq[i])
            index.append(i)
    return np.asarray(periods)


def get_frequency_index(freq, periods):
    # Return frequency index corresponding to the input periods
    # Input:
    # freq: Frequency array where index have to be extracted from
    # periods: Periods of analysis
    # Output:
    # index: List of index of freq array corresponding to periods
    index = [0 for i in range(len(periods))]

    count_periods = 0
    for i in range(int(np.floor(len(freq) / 2))):
        if freq[i + 1] > periods[count_periods]:
            index[count_periods] = i
            count_periods += 1
            if count_periods == len(periods):
                break
    if count_periods < len(periods):
        index = index[:count_periods]

    return index

def nextpow2(i):
    # Getting next power of 2, greater than i
    # Input:
    # i: number you wish to get the next power of 2
    # outputs:
    # n: power of 2 above i
    n = 1
    while n < i:
        n *= 2
    return n

def make_pow2(time, data):
    # Make signal length a power of 2 for FFT purposes
    # Input:
    # time: time array
    # data: signal array
    # Outputs:
    # time_2: time array padded with zeros so its length its a power of 2
    # data_2: data array padded with zeros so its length its a power of 2
    nfft = nextpow2(len(data))
    N = nfft - len(data)
    if not N == 0:
        data_2 = data.tolist()
        time_2 = time.tolist()
        data_2.extend([0 for i in range(N)])
        time_2.extend([time[-1] + (time[1] - time[0]) * i for i in range(N)])

    return np.asarray(time_2), np.asarray(data_2)

def get_power_spectra(sample_freq, data):
    # Get robust power spectra based on multitaper method
    # Input:
    # sample_freq: sample frequency of input data
    # data: Input data from which the power spectrum has to be derived
    # Output:
    # freq: frequency array of power spectrum
    # spec: robust power spectra
    nfft = 1024
    tbw = 3
    [tapers, eigen] = spectrum.dpss(nfft, tbw, 1)
    #res = spectrum.pmtm(data, e=tapers, v=eigen, show=False)
    amp, weights, _ = spectrum.pmtm(data, NW=3, show=False)
    freq = np.linspace(0, 1, len(amp[0])) * sample_freq

    spec = [0. for i in range(len(amp[0]))]

    for i in range(len(amp[0])):
        for j in range(len(amp)):
            spec[i] += amp[j][i] * weights[i][j]

    return freq, np.abs(spec)

def get_harmonics_amplitude(periods, freq, spectrum):
    # Get amplitude of harmonics
    # Input:
    # sample_freq: sample frequency of input data
    # data: Input data from which the power spectrum has to be derived
    # Output:
    # freq: frequency array of power spectrum
    # spec: robust power spectra

    df_log = 0.1 # width around harmonic for amplitude estimation
    df_log = 5 # for now, use indexes
    ps_values = [0. for i in range(len(periods))]

    for i in range(len(periods)):
        try:
            ps_values[i] = np.max([spectrum[int(periods[i]) - df_log:int(periods[i]) + df_log]])
        except ValueError:
            ps_values[i] = np.nan
    return ps_values

def get_distances(utm):
    # Calculates distances between two nodes
    # Input:
    # utm: list of 2 lists. utm[0]: list of easting; utm[1]: list of northing
    # Output:
    # distances: array of size(len(utm[0]), len(utm[0])) containing distances between nodes
    distances = np.asarray([[0. for i in range(len(utm[0]))] for j in range(len(utm[0]))])
    for i in range(len(utm[0])):
        distances[i, :] = np.asarray([np.sqrt((utm[0][i] - utm[0][j]) ** 2 + (utm[1][i] - utm[1][j]) ** 2) for j in range(len(utm[0]))])

    return distances

def get_power_ratio(values):
    # Calculates power ratio between two nodes
    # Input:
    # values: list of power for harmonics
    # Output:
    # ratios: array of size(len(values), len(values)) containing power ratio between nodes
    ratios = np.asarray([[0. for i in range(len(values))] for j in range(len(values))])
    for i in range(len(values)):
        ratios[i, :] = np.asarray([values[i] / values[j] for j in range(len(values))])

    return ratios

def remove_outliers(ps_values, periods):
    # Removes obvious outlier from harmonic selection based on variation in power and period
    # Input:
    # ps_values: values for power spectra
    # periods: frequencies corresponding to ps_values
    # Output:
    # ps_values_filt, periods_filt: list of filtered ps_values and periods

    error = [(ps_values[i + 1] - ps_values[i]) / ps_values[i] \
             for i in range(len(periods) - 1)]
    error_period = [(np.log10(periods[i + 2]) - np.log10(periods[i + 1])) / (np.log10(periods[i + 1]) - np.log10(periods[i])) \
             for i in range(len(periods) - 2)]

    ps_values_filt = []
    periods_filt = []
    for i in range(len(error)):
        if i < len(error) - 2:
            if error[i] < 5 and error_period[i] < 3:
                ps_values_filt.append(ps_values[i])
                periods_filt.append(periods[i + 1])
        else:
            if error[i] < 5:
                ps_values_filt.append(ps_values[i])
                periods_filt.append(periods[i])
    return ps_values_filt, periods_filt

def get_fft(time, data):
    # Get fft of 1D time series
    # Input:
    # time: array of time values
    # data: array of values for the 1D time series
    # Output:
    # freq: array of frequencies for spectral analysis
    # amp: Amplitude of power spectra
    # phase: phase of power spectra

    t_inc = time[1] - time[0]
    t_inc = 1. / 150.0
    if not np.mod(len(data), 2) == 0:
        nfft = nextpow2(len(data))
    else:
        nfft = len(data)
    spectrum = fftpack.fft(data, nfft)
    amp = np.abs(spectrum)
    phase = np.arctan2(spectrum.imag, spectrum.real)
    freq = np.linspace(0, 1, nfft) / t_inc

    return freq, amp, phase
