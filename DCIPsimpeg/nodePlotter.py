import matplotlib.pyplot as plt
import numpy as np
import DCIPtools as DCIP
from multiprocessing import Pool

# ============================ User area ==============================
# ver 0.1
# set variables
fname = "C:/Users/johnk/Projects/Karst/Lusha_North-QC.DAT"
# load the data file
patch = DCIP.loadDias(fname)

for rdg in range(len(patch.readings)):
    patch.readings[rdg].createNodeDB()
    for idx in range(len(patch.readings[rdg].node_db)):
        x = patch.readings[rdg].node_db[idx][0]
        y = patch.readings[rdg].node_db[idx][1]
        plt.plot(x, y, '*')
        plt.text(x, y, '%s' % patch.readings[rdg].node_db[idx])
    plt.show()


