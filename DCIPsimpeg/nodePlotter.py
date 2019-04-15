import matplotlib.pyplot as plt
import numpy as np
import DCIPtools as DCIP
from multiprocessing import Pool

# ============================ User area ==============================
# ver 0.1
# set variables
fname = "E:/Projects/debug/Mark/MMG_100m.DAT"
outname = "E:/Projects/debug/Mark/MMG_100m.log"
# load the data file
patch = DCIP.loadDias(fname)
patch.checkContinuityInNodeFiles(path=outname)


# print(nodes)
# rdg = 0;
# patch.readings[rdg].createNodeDB()
# patch.readings[rdg + 1].createNodeDB()
# print(patch.readings[rdg].node_ids)
# print(patch.readings[rdg + 1].node_ids)
# for idx in range(len(patch.readings[rdg].node_db)):
# for idx in range(len(patch.readings[rdg].node_db)):
#     x = patch.readings[rdg].node_locs[idx][0]
#     y = patch.readings[rdg].node_locs[idx][1]
#     plt.plot(x, y, '*')
#     plt.text(x, y, '%s' % patch.readings[rdg].node_db[idx])
# plt.show()


