import matplotlib.pyplot as plt
import numpy as np
import DCIPtools as DCIP
from multiprocessing import Pool

if __name__ == '__main__':
    # ============================ User area ==============================
    # ver 0.2
    # set variables
    fname = "E:/Projects/debug/Gunjan/All_50-400_CC.DAT"
    outname = "E:/Projects/debug/Gunjan/All_50-400_CC_decimate.DAT"
    # load the data file
    patch = DCIP.loadDias(fname)
    cnt = 0
    cnt_tot = 0
    # deg = []
    # diffs = []
    # bin1 = []
    # bin2 = []
    # bin3 = []
    # bin4 = []
    min_dp_size = patch.getSmallestDipoleSize()
    print("min dipole size: {0}".format(min_dp_size))
    for rdg in range(len(patch.readings)):
        for dp in range(len(patch.readings[rdg].Vdp)):
            cnt_tot+=1
            x0 = patch.readings[rdg].Vdp[dp].Rx1East
            y0 = patch.readings[rdg].Vdp[dp].Rx1North
            x1 = patch.readings[rdg].Vdp[dp].Rx2East
            y1 = patch.readings[rdg].Vdp[dp].Rx2North
            x2 = patch.readings[rdg].Idp.Tx1East
            y2 = patch.readings[rdg].Idp.Tx1North

            dp_size = ((x1 - x0)**2 + (y1 - y0)**2)**0.5
            tx_rx = ((x2 - x0)**2 + (y2 - y0)**2)**0.5
            separation = np.asarray([x0 - x1, y0 - y1])
            rx1 = [x0, y0]
            tx_rx = ((x2 - rx1[0])**2 + (y2 - rx1[1])**2)**0.5
            midpoint_rx = rx1 + separation / 2.
            if (separation[1] == 0):
                phi = 90
            else:
                phi = np.arctan(separation[0] / separation[1])
                phi = phi * 180 / np.pi
            tx_rx = ((x2 - midpoint_rx[0])**2 + (y2 - midpoint_rx[1])**2)**0.5
            if 45 > phi > -45:
                if 1200 > tx_rx > 500:
                    if dp_size > (7 * min_dp_size):
                        cnt+=1
                        # print(eta * 180 / np.pi)
                        # plt.plot([patch.readings[rdg].Vdp[dp].Rx1East,
                        #      patch.readings[rdg].Vdp[dp].Rx2East],
                        #      [patch.readings[rdg].Vdp[dp].Rx1North,
                        #      patch.readings[rdg].Vdp[dp].Rx2North], '-o')
                        # plt.plot(patch.readings[rdg].Idp.Tx1East,
                        #      patch.readings[rdg].Idp.Tx1North, 'dk')
                        # bin4.append(dp_size)
                    else:
                        patch.readings[rdg].Vdp[dp].flagRho = "Reject"
                        patch.readings[rdg].Vdp[dp].flagMx = "Reject"
                elif 500 > tx_rx > 200:
                    if (7 * min_dp_size) > dp_size > (4 * min_dp_size):
                        cnt+=1
                        # print(eta * 180 / np.pi)
                        # plt.plot([patch.readings[rdg].Vdp[dp].Rx1East,
                        #      patch.readings[rdg].Vdp[dp].Rx2East],
                        #      [patch.readings[rdg].Vdp[dp].Rx1North,
                        #      patch.readings[rdg].Vdp[dp].Rx2North], '-o')
                        # plt.plot(patch.readings[rdg].Idp.Tx1East,
                        #      patch.readings[rdg].Idp.Tx1North, 'dk')
                    else:
                        patch.readings[rdg].Vdp[dp].flagRho = "Reject"
                        patch.readings[rdg].Vdp[dp].flagMx = "Reject"
                        # bin3.append(dp_size)
                elif 200 > tx_rx > 0:
                    if (4 * min_dp_size) > dp_size:
                        cnt+=1
                        # print(eta * 180 / np.pi)
                        # plt.plot([patch.readings[rdg].Vdp[dp].Rx1East,
                        #      patch.readings[rdg].Vdp[dp].Rx2East],
                        #      [patch.readings[rdg].Vdp[dp].Rx1North,
                        #      patch.readings[rdg].Vdp[dp].Rx2North], '-o')
                        # plt.plot(patch.readings[rdg].Idp.Tx1East,
                        #      patch.readings[rdg].Idp.Tx1North, 'dk')
                    else:
                        patch.readings[rdg].Vdp[dp].flagRho = "Reject"
                        patch.readings[rdg].Vdp[dp].flagMx = "Reject"
                        # bin2.append(dp_size)
            else:
                cnt+=1
                # print(eta * 180 / np.pi)
                # plt.plot([patch.readings[rdg].Vdp[dp].Rx1East,
                #      patch.readings[rdg].Vdp[dp].Rx2East],
                #      [patch.readings[rdg].Vdp[dp].Rx1North,
                #      patch.readings[rdg].Vdp[dp].Rx2North], '-o')
                # plt.plot(patch.readings[rdg].Idp.Tx1East,
                #      patch.readings[rdg].Idp.Tx1North, 'dk')
                # patch.readings[rdg].Vdp[dp].flagRho = "Reject"
                # patch.readings[rdg].Vdp[dp].flagMx = "Reject"

    print(cnt, cnt_tot)
    # plt.ylim([patch.readings[rdg].Idp.Tx1North - 50, patch.readings[rdg].Idp.Tx1North + 50])
    # plt.axes().set_aspect('equal', 'datalim')
    # plt.show()
    strt_time = patch.window_start
    stop_time = patch.window_end
    patch.writeColeColeSEDat(outname, strt_time[0],
                             stop_time[stop_time.size - 1])
    # plt.plot(deg, '*')
    # plt.show()
