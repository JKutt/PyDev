import matplotlib.pyplot as plt
import numpy as np
import DCIPtools as DCIP

# ============================ User area ==============================
# set variables
fname = "/Users/juan/Documents/testData/B2_15_125-Loc.DATT"
outname = "/Users/juan/Documents/testData/B2_15_125-Loc-mod.DAT"
# load the data file
patch = DCIP.loadDias(fname)
G = []
for rdg in range(len(patch.readings)):
    for dp in range(len(patch.readings[rdg].Vdp)):
        G.append(patch.reading[rdg].Vdp[dp].getGeometricFactor(
                 patch.reading[rdg].Idp))

plt.hist(G)
plt.show()
