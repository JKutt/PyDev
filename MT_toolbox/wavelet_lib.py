import numpy as np
import spectrum
import pylab as plt
import scipy.linalg as slin
import time
from scipy.stats import rayleigh, beta
import scipy
import spectrum_lib

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


def local_maxima(wavelet_array):

    """
    # Automatic detection of local maxima
    # Input:
    # wavelet_array: modulus of wavelet analysis
    # Output:
    # location:
    """

    location = np.full(wavelet_array.shape, False)
    location_2 = np.full(wavelet_array.shape, 1)
    bool_1 = wavelet_array[:, 1:] > wavelet_array[:, :-1]
    bool_2 = wavelet_array[:, :-1] > wavelet_array[:, 1:]
    location[:, 1:-1] = bool_1[:, :-1] & bool_2[:, 1:]

    return location


def dyadic_frequencies(f_0, f_1, n_int):

    """
    # Wavelet analysis using Morlet
    # Input:
    # data: 1D time series for wavelet analysis
    # freq: Sample freq for data time series
    # freq_analysis: array of frequencies to analyze (will be turned into wavelet scaling)
    # Output:
    # wavelet_array: 2D array for wavelet analysis
    """

    f_0 = infpow2(f_0)
    f_1 = nextpow2(f_1)

    freq = []
    for i in range(getpow2(f_0, f_1)):
        f_init = f_0 * 2 ** i
        freq.extend([f_init * 2 ** (j / n_int) for j in range(n_int)])

    return freq

def get_scales(frequencies, w0):

    scales = [(w0 + np.sqrt(2 + w0 ** 2)) / 4 / np.pi / f_ana * np.pi * 2. for f_ana in frequencies]
    return scales

def wavelet_analysis(data, freq, freq_analysis, central_pulsation):

    """
    # Wavelet analysis using Morlet
    # Input:
    # data: 1D time series for wavelet analysis
    # freq: Sample freq for data time series
    # freq_analysis: array of frequencies to analyze (will be turned into wavelet scaling)
    # Output:
    # wavelet_array: 2D array for wavelet analysis
    """
    s_time = time.time()
    # Morlet parameters
    w0 = central_pulsation
    nfft = nextpow2(len(data))

    # Having power of 2 length for numerical speed (fftpack)
    data_fft = np.fft.fft(data, nfft)#.reshape(len(data))
    freq_array = np.linspace(0, 1, nfft) * freq

    # scales computation (morphing frequencies into Morlet scale)
    scales = get_scales(freq_analysis, w0)

    # More elegant wavelet cooking recipe
    wavelet_array = np.asarray([(2 * np.pi * scale / freq) ** (1. / 2.) * np.pi ** (- 1. / 4.) * np.exp(-1 * ((scale * freq_array - w0) ** 2) / 2) * np.complex(1, 0)
                     for scale in scales])

    # Removing negative frequencies to get the complex morlet wavelet
    wavelet_array[:, int(nfft // 2) + 1:] = np.complex(0, 0)
    # Zero mean for admissibility
    wavelet_array[:, 0] = np.complex(0, 0)
    # Correlation with signal
    wavelet_array *= np.tile(data_fft, (len(scales), 1))
    # Getting IFFT
    wavelet_array = np.abs(np.fft.ifft(wavelet_array, len(data), axis=1))

    """
    # slow loopy way
    # Initializing wavelet array
    wavelet_array = np.zeros((len(data), len(scales)))
    for scale, ind in zip(scales, range(len(scales))):
        morlet = (2 * np.pi * scale / freq) ** (1. / 2.) * np.pi ** (- 1. / 4.) * np.exp(-1 * ((scale * freq_array - w0) ** 2) / 2) * np.complex(1, 0)

        morlet[int(nfft // 2) + 1:] = [0. for i in range(len(morlet[int(nfft // 2) + 1:]))]
        morlet[0] = 0
        wavelet_array[:, ind] = np.abs(np.fft.ifft(morlet * data_fft, len(data)))
        #print(np.fft.ifft(morlet * data_fft, len(data)))
    """

    print("\t[INFO] Elapsed time for wavelet calculation: " + str(time.time() - s_time) + ' seconds.')

    return wavelet_array


def morlet_kernel(central_pulsation, scale, scales, sample_freq):

    """
    # Return the morlet correlation kernel (Maraun, 2006)
    # Source: What can we learn from climate data? Methods for fluctuation, Time/Scale and phase analysis, PhD thesis
    # Input:
    # central_pulsation: central pulsation (w0 in equations) of Morlet wavelet
    # scale: center scale for Morlet kernel
    # scales: other set of scales to compute the correlation kernel
    # sample_freq: sampling frequency
    # Output:
    # kernel: normalized morlet kernel (absolute value)
    """

    time = np.linspace(1, 32768, 32768) / sample_freq # Creating time array
    t0 = time[int(len(time) / 2)]     # Putting the kernel at the center

    # Computing kernel, scale by scale
    kernel = np.zeros((len(scales), len(time)))
    for i in range(len(scales)):
        kernel[i, :] = np.abs(np.sqrt((2. * scales[i] * scale) / (scales[i] ** 2 + scale ** 2)) *\
                        np.exp(np.complex(0, 1) * central_pulsation * (scales[i] + scale) / (scales[i] ** 2 + scale ** 2) * (t0 - time)) *\
                        np.exp(-0.5 * ((t0 - time) ** 2 + central_pulsation ** 2 * (scales[i] - scale) ** 2) / (scales[i] ** 2 + scale ** 2)))

    #plt.imshow(kernel / np.max(np.max(kernel)), aspect='auto')
    #plt.show()
    return kernel / np.max(np.max(kernel))


def correlation_distance(kernel, criterion):

    """
    # Return the correlation distance defined by input criteria (the bigger the criteria, the higher the correlation is)
    # Input:
    # kernel: morlet kernel
    # criterion: 0 (no correlation)< alpha < 1 (total correlation)
    # Output:
    # distance:
    """

    kernel[kernel < criterion] = 0 # Removing low correlation kernel for distance computation
    length = np.max([len(kernel[i, :][kernel[i, :] > 0]) for i in range(len(kernel))])

    return int(length)


def morlet_correlation_distance(central_pulsation, frequencies, sample_freq, criterion):

    scales = [(central_pulsation + np.sqrt(2 + central_pulsation ** 2)) / 4 / np.pi / f_ana * np.pi * 2. for f_ana in frequencies]
    distances = [0 for scale in scales]
    for (scale, i) in zip(scales, range(len(distances))):
        kernel = morlet_kernel(central_pulsation, scale, scales, sample_freq)
        distances[i] = correlation_distance(kernel, criterion)

    return distances


def chains_maxima(maxima, distances, frequencies, limits):

    #distances = [10 for i in range(len(distances))]
    n_scales = len(maxima[:, 0])
    n_time = len(maxima[0, :])
    chains = []
    start_frequency = 2 ** ((np.log2(limits[1]) + np.log2(limits[0])) / 2)
    tmp = np.abs(start_frequency - np.asarray(frequencies))
    ind = [i for i in range(len(frequencies)) if np.min(tmp) == tmp[i]][0]
    tmp = np.abs(limits[0] - np.asarray(frequencies))
    ind_min = [i for i in range(len(frequencies)) if np.min(tmp) == tmp[i]][0]
    tmp = np.abs(limits[1] - np.asarray(frequencies))
    ind_max = [i for i in range(len(frequencies)) if np.min(tmp) == tmp[i]][0]

    potential_chains = [i for i in range(n_time) if maxima[ind, i]]
    for init in potential_chains:
        tmp_chain = [init]                  # Initializing at mid-frequency
        for i in range(ind + 1, len(frequencies)):         # chaining up
            chained = 0                     # Criteria to complete chains with nan if chaining failed
            # exploring on the time axis, further and further away every iteration
            for j in range(distances[i]):
                if init + j < n_time:
                    if maxima[i, init + j]: # If close local maxima
                        chained = 1
                        tmp_chain.append(init + j)
                        break               # Moving on to next scale
                if init - j >= 0:
                    if maxima[i, init - j]: # If close local maxima
                        chained = 1
                        tmp_chain.append(init - j)
                        break               # Moving on to next scale
            if j == distances[i] - 1 and chained == 0: # Completing with nan if chaining failed
                tmp_chain.extend([np.nan for k in range(i, len(frequencies))])
                break

        for i in list(reversed(range(0, ind))):            # chaining down
            chained = 0                     # Criteria to complete chains with nan if chaining failed
            # exploring on the time axis, further and further away every iteration
            for j in range(distances[i]):
                if init + j < n_time:
                    if maxima[i, init + j]: # If close local maxima
                        chained = 1
                        tmp_chain.insert(0, init + j)
                        break               # Moving on to next scale
                if init - j >= 0:
                    if maxima[i, init - j]: # If close local maxima
                        chained = 1
                        tmp_chain.insert(0, init - j)
                        break               # Moving on to next scale
            if j == distances[i] - 1 and chained == 0: # Completing with nan if chaining failed
                for k in list(reversed(range(0, i + 1))):
                    tmp_chain.insert(0, np.nan)
                break

        if not np.isnan(tmp_chain[ind_max]) and not np.isnan(tmp_chain[ind_min]):
            chains.append(tmp_chain)
    return chains

def get_events(output, input, ref, distances):

    events_input = []
    time_input = [[], []]
    time_ref = [[], []]
    events_ref = []
    chains_out = []

    for i in range(len(input)):
        for j in range(len(input[i])):
            time_input[i].append(int(np.median([input[i][j][k] for k in range(len(distances)) if not np.isnan(input[i][j][k])])))

    for time in time_input[0]:
        diff = np.abs(time - np.asarray(time_input[1]))
        if np.min(diff) <= 10:
            events_input.append(time)

    for i in range(len(ref)):
        for j in range(len(ref[i])):
            time_ref[i].append(int(np.median([ref[i][j][k] for k in range(len(distances)) if not np.isnan(ref[i][j][k])])))

    for time in time_ref[0]:
        diff = np.abs(time - np.asarray(time_ref[1]))
        if np.min(diff) <= 10:
            events_ref.append(time)


    events_mag = []
    # Comparison between local and remote mag
    time_dif = [[np.abs(events_input[i] - events_ref[j])
                    for i in range(len(events_input))]
                    for j in range(len(events_ref))]

    for i in range(len(events_input)):
        for j in range(len(events_ref)):
            if time_dif[j][i] <= np.max([np.min(distances), 5]):
                events_mag.append(events_input[i])

    events = [[] for i in range(len(output))]
    for i in range(len(output)):
        time_output = []
        # comparing chains between input channels
        for j in range(len(output[i])):
            time_output.append(int(np.median([output[i][j][k] for k in range(len(distances)) if not np.isnan(output[i][j][k])])))

        for time in time_output:
            diff = np.abs(time - np.asarray(events_mag))
            if np.min(diff) <= np.max([np.min(distances), 5]):
                events[i].append(time)

        # comparing chains with output and input
    return [events_mag, events_mag]


def global_wavelet_spectrum(coefficients):

    nf = len(coefficients[:, 0])
    nt = len(coefficients[0, :])

    global_spectrum = np.asarray([0. for i in range(nf)])

    for i in range(nf):
        tmp = []
        tmp = coefficients[i, :]
        tmp = tmp[~np.isnan(tmp)]
        global_spectrum[i] = np.sum(np.abs(tmp)**2)
        global_spectrum[i] = 1. / len(coefficients[i, :]) * global_spectrum[i]

    return global_spectrum


def significant_coefficients(signal, wavelets, sample_freq, central_pulsation, freq_analysis, level):

    surrogate = spectrum_lib.fourier_surrogate(signal, sample_freq)
    wavelets_surrogate = wavelet_analysis(surrogate, sample_freq, freq_analysis, central_pulsation)
    global_spectrum_surrogate = global_wavelet_spectrum(wavelets_surrogate)

    if level == 0.7:
        chi = 2.408
    if level == 0.8:
        chi = 3.219
    if level == 0.9:
        chi = 4.605
    if level == 0.95:
        chi = 5.991
    if level == 0.975:
        chi = 7.378
    if level == 0.99:
        chi = 9.210
    if level == 0.995:
        chi = 10.59

    for j in range(len(wavelets[:, 0])):
        wavelets[j, np.abs(wavelets[j, :]) ** 2 < global_spectrum_surrogate[j] * 0.5 * chi] = np.nan

    return wavelets
