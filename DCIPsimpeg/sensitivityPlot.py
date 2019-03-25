import numpy as np
import DCIPtools as DCIP
import matplotlib.pyplot as plt
from multiprocessing import Pool


def fastJknockout(reading):
    sw_s = np.asarray([389302, 3278960])
    ne_s = np.asarray([389368, 3279100])
    dx = 12.5
    J_ = np.zeros((91, 161, 19))
    current = reading.Idp

    for dp in range(len(reading.Vdp)):
        J_ += reading.Vdp[dp].calcSensitivity(current, ne, sw, dx, sw_search=sw_s, ne_search=ne_s)
    
    return reading

# ============================ User area ==============================
# ver 0.1
# set variables
fname = "C:\\Users\\johnk\\Projects\\Karst\\Lusha_North-QC.DAT"
outname = "C:\\Users\\johnk\\Projects\\Karst\\Lusha_North-QC-Jflag.DAT"
# load the data file
patch = DCIP.loadDias(fname)

rx_locs = patch.getDipoles()
tx_locs = patch.getSources()
print("min E: {0} & max E: {1}".format(np.min(rx_locs[:, 0]),
      np.max(rx_locs[:, 0])))
print("min N: {0} & max N: {1}".format(np.min(rx_locs[:, 1]),
      np.max(rx_locs[:, 1])))
print("min Elev: {0} & max Elev: {1}".format(np.min(rx_locs[:, 2]),
      np.max(rx_locs[:, 2])))
sw = np.array([np.min(rx_locs[:, 0]), np.min(tx_locs[:, 1]),
              np.min(rx_locs[:, 2]) - 200])
ne = np.array([np.max(rx_locs[:, 0]), np.max(rx_locs[:, 1]),
              np.max(rx_locs[:, 2])])
# rdg = 10
# dp = 15
sw_s = np.asarray([388712, 3279214])
ne_s = np.asarray([388816, 3279319])
# sw_s = np.asarray([389302, 3278960])
# ne_s = np.asarray([389368, 3279100])
dx = 12.5
# pool = Pool()
# # run the inversion scheme
# pro = pool.map(fastJknockout, patch.readings)
# patch.readings = pro

J_ = np.zeros((91, 161, 19))

for rdg in range(len(patch.readings)):
	current = patch.readings[rdg].Idp
	for dp in range(len(patch.readings[rdg].Vdp)):
		J_ += patch.readings[rdg].Vdp[dp].calcSensitivity(current,
		                                                  ne, sw, dx,
		                                                  sw_search=sw_s,
		                                                  ne_search=ne_s)
# current = patch.readings[rdg].Idp
# J_ += patch.readings[rdg].Vdp[dp].calcSensitivity(current,
# 		                                          ne, sw, dx,
# 		                                          sw_search=sw_s,
# 		                                          ne_search=ne_s)
start_time = 40
end_time = 2000
start_inds = (patch.window_start >= start_time)
stop_inds = (patch.window_end <= end_time)
strt_time = patch.window_start[start_inds]
stop_time = patch.window_end[stop_inds]
# write the data to file
print(strt_time)
patch.writeColeColeSEDat(outname, strt_time[0],
                         stop_time[stop_time.size - 1])
# X = np.arange(sw[0], ne[0], dx)
# Y = np.arange(sw[1], ne[1], dx)
# x, y = np.meshgrid(X, Y)
# z = J_[:, :, 3]
# # print(np.min(x[a]), np.max(x[a]))
# # print(np.min(y[a]), np.max(y[a]))
# plt.contourf(x, y, z, cmap='viridis')
# plt.colorbar()
# tx1 = [patch.readings[rdg].Idp.Tx1East, patch.readings[rdg].Idp.Tx1North]
# rx1 = [patch.readings[rdg].Vdp[dp].Rx1East, patch.readings[rdg].Vdp[dp].Rx1North]
# rx2 = [patch.readings[rdg].Vdp[dp].Rx2East, patch.readings[rdg].Vdp[dp].Rx2North]
# plt.plot(tx1[0], tx1[1], 'r*')
# # plt.plot(sbx, sby, 'k*')
# plt.plot([sw_s[0], ne_s[0]], [sw_s[1], ne_s[1]], 'kd')
# plt.plot(rx1[0], rx1[1], 'bo')
# plt.plot(rx2[0], rx2[1], 'go')
# plt.axes().set_aspect('equal', 'datalim')
# plt.show()
