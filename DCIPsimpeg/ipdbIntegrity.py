import matplotlib.pyplot as plt
import numpy as np
import DCIPtools as DCIP
# from multiprocessing import Pool

# ============================ User area ==============================
# ver 0.1
# set variables
fname = "E:/Projects/China/Karst/IPDB/Lusha_North-QC.DAT"
outname = "E:/Projects/China/Karst/IPDB/Lusha_North-QC.log"
# utm_csv = 
# load the data file
patch = DCIP.loadDias(fname)

# # 1) check the continuity of the files used
# """
#     goes through database and creates node database for each reading
#     and checks that file number id's are consecutive

#     OUTPUT: log file indicating discontinuities
# """
# patch.checkContinuityInNodeFiles(path=outname)


# # 2) checks records in log file to records actually processed 
# """ 
#     runs through database and compares availble records to
#     records log
# """
# rec_log_fpath = "E:/Projects/China/Karst/IPDB/B2/Log/Recordings_0.txt"
# patch.checkForMissingRecords(record_log=rec_log_fpath)


# # 3) plot node locations + Node ID for a defined reading
# """
#     runs through each reading and plots the avaible nodes and dipole
#     polarity
# """
# for rdg in range(len(patch.readings)):
#     patch.readings[rdg].createNodeDB()
#     for idx in range(len(patch.readings[rdg].node_db)):
#         x = patch.readings[rdg].node_locs[idx][0]
#         y = patch.readings[rdg].node_locs[idx][1]
#         plt.plot(x, y, 'k*')
#         plt.text(x, y, '%s' % patch.readings[rdg].node_db[idx])
#     plt.title("rec: {0}".format(patch.readings[rdg].MemNumber))
#     for dp in range(len(patch.readings[rdg].Vdp)):
#         rx1x = patch.readings[rdg].Vdp[dp].Rx1East
#         rx1y = patch.readings[rdg].Vdp[dp].Rx1North
#         rx2x = patch.readings[rdg].Vdp[dp].Rx2East
#         rx2y = patch.readings[rdg].Vdp[dp].Rx2North
#         In = patch.readings[rdg].Idp
#         if patch.readings[rdg].Vdp[dp].calcRho(In) > 0:
#             plt.plot([rx1x, rx2x], [rx1y, rx2y], '-ob')
#         else:
#             plt.plot([rx1x, rx2x], [rx1y, rx2y], '-or')
#     plt.show()

# 4) histograms of vital data
"""
   runs through database and creates histogram of the data
   vitals such as rho, Mx, Vp
"""
fig = plt.figure()
ax1 = fig.add_subplot(241)
ax2 = fig.add_subplot(242)
ax3 = fig.add_subplot(243)
ax4 = fig.add_subplot(244)
ax5 = fig.add_subplot(245)
ax6 = fig.add_subplot(246)
ax7 = fig.add_subplot(247)
ax8 = fig.add_subplot(248)

ax1.set_title(r'$\rho$')
ax2.set_title("error V/I")
ax3.set_title("K")
ax4.set_title("V")
ax5.set_title(r'$\eta$')
ax6.set_title("Null")
ax7.set_title("I")
ax8.set_title("a")

rho = patch.getApparentResistivity()
vp_i_error = patch.getVpDivInErrors()
gk = patch.getGeometricFactor()
voltages = patch.getVoltages()
mx = patch.getApparentChageability()
null = patch.getNullFactors()
current = patch.getCurrents()
a = patch.getAspacings()

ax1.hist(rho, 50)
ax2.hist(vp_i_error, 10)
ax3.hist(gk)
ax4.hist(voltages, 50)
ax5.hist(mx,)
ax6.hist(null)
ax7.hist(current)
ax8.hist(a)

plt.show()
