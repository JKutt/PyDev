import numpy as np
import spectrum
import pylab as plt
import scipy.linalg as slin
import time
from scipy.stats import rayleigh, beta
import scipy

def nextpow2(length):

    i = 0
    while 2 ** i < length:
        i += 1

    return 2 ** i


def infpow2(length):

    i = 0
    while 2 ** i <= length:
        i += 1

    return 2 ** (i - 1)


def getpow2(length_1, length_2):

    i = 0
    while length_1 * 2 ** i <= length_2:
        i += 1

    return i - 1

def dyadic_frequencies(f_0, f_1, n_int):

    f_0 = infpow2(f_0)
    f_1 = nextpow2(f_1)

    freq = []
    for i in range(getpow2(f_0, f_1)):
        f_init = f_0 * 2 ** i
        freq.extend([f_init * 2 ** (j / n_int) for j in range(n_int)])

    return freq


def wavelet_analysis(data, freq, freq_analysis):

    # Morlet parameters
    w0 = 5
    nfft = nextpow2(len(data))

    # Having power of 2 length for numerical speed
    data_fft = np.fft.fft(data, nfft)
    freq_array = np.linspace(0, 1, nfft) * freq


    # scales computation
    scales = [(w0 + np.sqrt(2 + w0 ** 2)) / 4 / np.pi / f_ana * np.pi * 2. for f_ana in freq_analysis]
    wavelet_analysis = np.zeros((len(data), len(scales)))

    for scale, ind in zip(scales, range(len(scales))):
        morlet = (2 * np.pi * scale / freq) ** (1. / 2.) * np.pi ** (- 1. / 4.) * np.exp(-1 * ((scale * freq_array - w0) ** 2) / 2) * np.complex(1, 0)
        morlet[int(nfft // 2) + 1:] = [0. for i in range(len(morlet[int(nfft // 2) + 1:]))]
        morlet[0] = 0
        wavelet_analysis[:, ind] = np.abs(np.fft.ifft(morlet * data_fft, len(data)))
        #print(np.fft.ifft(morlet * data_fft, len(data)))

    #print()

    return wavelet_analysis
