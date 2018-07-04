import matplotlib.pyplot as plt
# import numpy as np
import JDataObject as Jdata

# define the file required for import
fileName = "C:/Users/johnk/Projects/teck/6050/6050-2.DAT"

patch = Jdata.loadDias(fileName)   # Create the patch from data file

rdg = 0                           # Source to plot
dpnum = []
for i in range(len(patch.readings[rdg].Vdp)):
    dpnum.append(patch.readings[rdg].Vdp[i].dipole)
    plt.plot(patch.window_center,
             patch.readings[rdg].Vdp[i].Vs /
             (patch.readings[rdg].Vdp[i].Vp / 1000.), '-o')
plt.legend(dpnum)
plt.title("All Decay")
plt.xlabel("time (ms)")
plt.ylabel("Voltage (mV)")
plt.show()
