
import DCIPtools as DCIP
import matplotlib.pyplot as plt
import numpy as np

################################################################
# Read binary time-series ======================================


def getData():
    xyz = open("/Users/juan/Documents/testData/6050_R159_OY_OC.xyz")
    t = []
    xt = []
    samples_half_t = 600.0
    for line in xyz:
        x, y = line.split()
        t.append(float(x))
        xt.append(float(y))
    xyz.close()
    xt = np.asarray(xt)
    num_half_T = np.floor(xt.size / samples_half_t)
    trim_length = num_half_T * samples_half_t
    xt = xt[0:int(trim_length)]
    xt = np.asarray(xt)
    return xt, num_half_T


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
# end read =====================================================


# stack the data ===================================================
xt, num_half_T = getData()
timeFrom, timeTo = getTime()
start_vp = 50                           # start of Vp calculation (%)
end_vp = 90                             # end of Vp calculation (%)
window = DCIP.createHanningWindow(num_half_T)   # creates filter window
window2 = DCIP.createSlepianWindow(int(num_half_T), 0.3)
window3 = DCIP.createChebyshevWindow(int(num_half_T), 500)
window4 = DCIP.createKaiserWindow(int(num_half_T), 54)
# window3 = DCIP.createHanningWindow(7)   # creates filter window
# window2 = DCIP.createBruteStackWindow(int(num_half_T))
# # print(window2.size, window.size)
tHK = DCIP.filterKernal(filtershape=window)     # creates filter kernal
tHK2 = DCIP.filterKernal(filtershape=window2)     # creates filter kernal
print("half T: {0} window: {1} Kernel: {2}".format(num_half_T, window.size, tHK.kernel.size))
# # eHK = DCIP.ensembleKernal(filtershape=window3,
# #                           number_half_periods=num_half_T)
dkernal = DCIP.decayKernal(num_windows=np.asarray(timeTo).size,
                           window_starts=np.asarray(timeFrom),
                           window_ends=np.asarray(timeTo),
                           window_weight=501,
                           window_overlap=0.99,
                           output_type="Vs")  # creates decay kernal
stack = tHK * xt                               # stack data
stack2 = tHK2 * xt                               # stack data
# # ens = eHK * xt
# # plt.plot(stack)
# # plt.show()
decay = dkernal * (tHK * xt)         # calculates the decay
decay2 = dkernal * (tHK2 * xt)         # calculates the decay
# # end  ======================================================

# # cole-cole fitting =========================================
# Vp = DCIP.getPrimaryVoltage(start_vp,
#                             end_vp,
#                             stack)      # calculate the Vp
# staticmethod
# # # end =======================================================
# # plt.plot(dkernal.getWindowCenters(), win_std)
# # plt.show()
# plot results
# fig = plt.figure(figsize=(10, 8))
# ax1 = plt.subplot(311)
# ax1.plot(xt, 'r')
# ax1.set_xlabel("num samples")
# ax1.set_ylabel("voltage (mV)")
# ax1.set_title("Raw time-series")
# ax4 = plt.subplot(323)
# ax4.plot(tHK.getFilterKernal, 'og')
# ax4.set_xlabel("num taps")
# ax4.set_ylabel("amplitude")
# amp = DCIP.getFrequnceyResponse(tHK.getFilterKernal)
# amp2 = DCIP.getFrequnceyResponse(tHK2.getFilterKernal)
# freqs = np.arange(0, amp.size) * (150. / window.size)
# ax5 = plt.subplot(324)
# ax5.plot(freqs, amp, 'm')
# ax5.plot(freqs, amp2, 'k')
# ax5.set_xlabel("frequency (Hz)")
# ax5.set_ylabel("amplitude")
# ax2 = plt.subplot(326)
# ax2.plot(dkernal.getWindowCenters(), decay, '-ko')
# ax2.plot(dkernal.getWindowCenters(), decay2, '-ro')
# # ax2.plot(dkernal.getWindowCenters(), vs, '-o')
# ax2.set_xlabel("time (ms)")
# ax2.set_ylabel("Voltage (mV)")
# ax2.set_title("Secondary Voltage (decay)")
# ax3 = plt.subplot(325)
# ax3.plot(tHK * -xt)
# ax3.set_xlabel("num samples")
# ax3.set_ylabel("Voltage (mV)")
# ax3.set_title("Stack (decay)")
fig = plt.figure(figsize=(10, 8))
ax1 = plt.subplot(311)
ax1.plot(xt, 'r')
ax1.set_xlabel("num samples")
ax1.set_ylabel("voltage (mV)")
ax1.set_title("Raw time-series")
ax4 = plt.subplot(312)
ax4.plot(tHK.getFilterKernal, 'og')
ax4.set_xlabel("num taps")
ax4.set_ylabel("amplitude")
plt.show()
