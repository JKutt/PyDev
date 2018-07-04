import matplotlib.pyplot as plt
import numpy as np
import JDataObject as Jdata

# define the file required for import
fileName = "C:/Users/johnk/Downloads/GridB_800_GridLoc.DAT"

patch = Jdata.loadDias(fileName)   # Create the patch from data file

rdg = 0                           # Source to plot
# i = 4
dpnum = []

for rdg in range(len(patch.readings)):
    for i in range(len(patch.readings[rdg].Vdp)):
        rx2 = np.asarray([patch.readings[rdg].Vdp[i].Rx2East,
                         patch.readings[rdg].Vdp[i].Rx2North])
        rx1 = np.asarray([patch.readings[rdg].Vdp[i].Rx1East,
                         patch.readings[rdg].Vdp[i].Rx1North])
        Tx1 = np.asarray([patch.readings[rdg].Idp.Tx1East,
                         patch.readings[rdg].Idp.Tx1North])
        Tx2 = np.asarray([266476.000, 7318971.000])

        separation = rx2 - rx1
        midpoint_rx = rx1 + separation / 2.
        txrx = midpoint_rx - Tx1
        # print(separation)
        alpha = (180.0 / np.pi) * np.arctan(txrx[1] / txrx[0])
        beta = (180.0 / np.pi) * np.arctan(separation[1] / separation[0])
        eta = alpha - beta
        # print(alpha, beta, eta, patch.readings[rdg].Vdp[i].dipole)
        plt.plot([patch.readings[rdg].Vdp[i].Rx1East,
                 patch.readings[rdg].Vdp[i].Rx2East],
                 [patch.readings[rdg].Vdp[i].Rx1North,
                 patch.readings[rdg].Vdp[i].Rx2North], '-o')
        plt.plot(patch.readings[rdg].Idp.Tx1East,
                 patch.readings[rdg].Idp.Tx1North, 'dk')
    # plt.plot([midpoint_rx[0], Tx1[0]], [midpoint_rx[1], Tx1[1]], 'd-r')
    # plt.legend(dpnum)
plt.axes().set_aspect('equal', 'datalim')
plt.title("All Dipoles")
plt.xlabel("Easting (m)")
plt.ylabel("Northing (m)")
plt.show()

