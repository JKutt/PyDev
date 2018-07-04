import matplotlib.pyplot as plt
import numpy as np
import JDataObject as Jdata

# define the file required for import
fileName = "/home/juan/Documents/DevProjects/testData/L5150-final.DAT"

patch = Jdata.loadDias(fileName)   # Create the patch from data file
rdg = 0
i = 15
rho = []
mx = []

# plt.plot(patch.window_center,
#          patch.readings[rdg].Vdp[i].Vs, '-o')

win_length = 7
side = (win_length - 1) / 2.
zero_pad = (win_length - 1) / 2.0
sides = zero_pad
kappa = 1.4826
eta = 1.0
win = np.zeros(win_length)
print(patch.readings[rdg].Vdp[i].Vs)
for index in range(len(patch.readings)):
    for dp in range(len(patch.readings[index].Vdp)):
        if 80 > patch.readings[index].Vdp[dp].Mx > 0:
            rho.append(patch.readings[index].Vdp[dp].Rho)
            mx.append(patch.readings[index].Vdp[dp].Mx)

local_mx = []
rdg = 0
dp = 4
for index in range(len(patch.readings)):
    for dp in range(len(patch.readings[index].Vdp)):
        
        dist = np.sqrt()

# for index in range(patch.window_center.size):
    # if (index - side) >= 0 and (index + side) <= (patch.window_center.size - 1):
    #     fr = int(index - side)
    #     to = int(index + side + 1)
    #     win = patch.readings[rdg].Vdp[i].Vs[fr:to]
    #     local_median = np.median(win)
    #     mean_abs_dev = np.median(np.abs(win - local_median))
    #     sigma = eta * mean_abs_dev
    #     lhs = np.abs(patch.readings[rdg].Vdp[i].Vs[i] - local_median)
    #     rhs = np.abs(eta * sigma)
    #     # print(lhs, rhs)
    #     if lhs > rhs:
    #         print("found outlier ", index)
    #         patch.readings[rdg].Vdp[i].Vs[index] = local_median

# plt.plot(patch.window_center,
#          patch.readings[rdg].Vdp[i].Vs, 'r-o')
print(np.median(mx), np.mean(mx))
plt.hist(mx)
# plt.legend(dpnum)
plt.title("Mx")
# plt.xlabel("time (ms)")
# plt.ylabel("Voltage (mV/V)")
plt.show()
