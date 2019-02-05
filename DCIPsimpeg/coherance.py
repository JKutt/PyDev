import DCIPtools as DCIP
import matplotlib.pyplot as plt
import numpy as np
from scipy import fftpack

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
xt2 = -xt  # np.random.rand(xt.size)
xt = DCIP.padNextPower2(xt)
xt2 = DCIP.padNextPower2(xt2)
freq, coh = DCIP.getCoherence(xt, xt2, 150)
x_corr = DCIP.getCrossCorrelation(xt, xt2)

# plt.plot(xt, 'r')
# plt.plot(xt2, 'b')
plt.plot(x_corr)
plt.show()
