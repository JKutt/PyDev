#!C:\\Users\\HugoLarnier\\Anaconda3\\python

from datfiles_lib import *
import csv
import matplotlib.pylab as plt
import matplotlib.lines as mlines
from glob import glob
import numpy as np
import DCIPtools as DCIP
import signalAnalysis as signal

def getTime():
    timeFrom = [2040., 2060., 2080., 2120., 2160., 2200.,
                2240., 2320., 2400.,
                2480., 2560., 2640.,
                2720., 2800., 2960.,
                3120., 3280., 3440.,
                3600., 3760.]
    timeTo = [2060., 2080., 2120., 2160., 2200., 2240.,
              2320., 2400., 2480., 2560., 2640., 2720.,
              2800., 2960., 3120., 3280., 3440.,
              3600., 3760., 3920.]
    return timeFrom, timeTo
timeFrom, timeTo = getTime()

#node = 'C:\\Users\\HugoLarnier\\Documents\\projects\\mmg_mcArthur\\L34\\L34\\DATA\\TXR\\D6000129.DAT'
node = 'C:\\Users\\HugoLarnier\\Documents\\projects\\mmg_mcArthur\\L41\\L41\\DATA\\E3\\E3001152.DAT'
node = 'C:\\Users\\HugoLarnier\\Documents\\projects\\Orano\\data\\SD0\\DATA\\AV\\AV000987.DAT'
fIn = open(node, 'r')
linesFIn = fIn.readlines()
fIn.close()

time, data = read_data(linesFIn)
data = data[1:]
time = time[1:]

sample_freq = 150.
sample_half_t = 600.0
time_s, data_s = signal.synthetic_current(sample_freq, 1)

periods = signal.get_maxima(1,
                            len(data_s),
                            np.linspace(0, 1, len(data_s)) * sample_freq,
                            np.abs(np.fft.fft(data_s)))
index_2 = signal.get_frequency_index(np.linspace(0, 1, len(data_s)) * sample_freq, periods)
ps_values_2 = signal.get_harmonics_amplitude(index_2,
                                            np.linspace(0, 1, len(data_s)) * sample_freq,
                                            np.abs(np.fft.fft(data_s)))

ps_values, periods = signal.remove_outliers(ps_values_2, periods)

fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, sharex=True)
ax1.plot(np.linspace(0, 1, len(data_s)) * sample_freq, np.abs(np.fft.fft(data_s)))
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.plot(periods, ps_values, 'ko')
ax2.plot(periods[:-1], [(ps_values[i + 1] - ps_values[i]) / ps_values[i] \
         for i in range(len(periods) - 1)], 'ro')
ax3.plot(periods, ps_values, 'ko')
ax4.plot(periods[:-2], [(np.log10(periods[i + 2]) - np.log10(periods[i + 1])) / (np.log10(periods[i + 1]) - np.log10(periods[i])) \
         for i in range(len(periods) - 2)], 'ro')
plt.show()

num_half_T = np.round(np.asarray(data).size / sample_half_t) - 2
trim_length = num_half_T * sample_half_t
data = np.asarray(data[0:int(trim_length)])
time = np.asarray(time[0:int(trim_length)])

freq2, amp2 = signal.get_power_spectra(sample_freq, data)

index = signal.get_frequency_index(freq2, periods)
ps_values = signal.get_harmonics_amplitude(index, freq2, amp2)

fig, (ax2, ax3) = plt.subplots(2, 1, sharex=True)
ax2.plot(periods[:-1], [periods[i + 1] - periods[i] for i in range(len(periods) - 1)], 'ro')
ax3.plot(periods, ps_values, 'ko')
plt.show()
fig, (ax1, ax2) = plt.subplots(2, 1)
ax1.plot(freq2, amp2, 'r')
ax1.plot(periods,
         ps_values, 'ok')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax2.plot(time, data)
plt.show()

start_vp = 50                           # start of Vp calculation (%)
end_vp = 90                             # end of Vp calculation (%)
window = DCIP.createHanningWindow(num_half_T)   # creates filter window
window3 = DCIP.createHanningWindow(7)   # creates filter window
window2 = DCIP.createBruteStackWindow(int(num_half_T))

tHK = DCIP.filterKernal(filtershape=window)     # creates filter kernal
tHK.kernel = tHK.kernel[1:-1]

# eHK = DCIP.ensembleKernal(filtershape=window3,
#                           number_half_periods=num_half_T)
dkernal = DCIP.decayKernal(num_windows=np.asarray(timeTo).size,
                           window_starts=np.asarray(timeFrom),
                           window_ends=np.asarray(timeTo),
                           window_weight=301,
                           window_overlap=0.95,
                           output_type="Vs")


stack = tHK * data                               # stack data
decay = dkernal * (tHK * np.asarray(data))         # calculates the decay
Vp = DCIP.getPrimaryVoltage(start_vp,
                            end_vp,
                            stack)      # calculate the Vp
print(Vp)
fig = plt.figure(figsize=(10, 8))
ax1 = plt.subplot(311)
ax1.plot(data, 'r')
ax1.set_xlabel("num samples")
ax1.set_ylabel("voltage (mV)")
ax1.set_title("Raw time-series")
ax4 = plt.subplot(323)
ax4.plot(tHK.getFilterKernal, 'og')
ax4.set_xlabel("num taps")
ax4.set_ylabel("amplitude")
amp = DCIP.getFrequencyResponse(window)
freqs = np.arange(0, amp.size) * (150. / window.size)
ax5 = plt.subplot(324)
ax5.plot(freqs, amp, 'm')
ax5.set_xlabel("frequency (Hz)")
ax5.set_ylabel("amplitude")
ax2 = plt.subplot(326)
ax2.plot(dkernal.getWindowCenters(), decay, '-ko')
# ax2.plot(dkernal.getWindowCenters(), decay, '.')
# ax2.plot(dkernal.getWindowCenters(), vs, '-o')
ax2.set_xlabel("time (ms)")
ax2.set_ylabel("Voltage (mV)")
ax2.set_title("Secondary Voltage (decay)")
ax3 = plt.subplot(325)
ax3.plot(tHK * -data)
ax3.set_xlabel("num samples")
ax3.set_ylabel("Voltage (mV)")
ax3.set_title("Stack (decay)")
plt.show()
