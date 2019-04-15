import matplotlib.pyplot as plt
import numpy as np
import DCIPtools as DCIP
import matplotlib.colors as colors

fname = "C:\\Users\\johnk\\Downloads\\HPX_Tintic_Swath1+2_AllDeliverables\\AllDeliverables\\Tintic_Swath2_Gradient_AllDipoles_v2\\Tintic_Swath2_Gradient_AllDipoles_v2.dat"
# load the data file
patch = DCIP.loadDias(fname)

for rec in range(len(patches.readings)):
rec = 5
print(patch.readings[rec].MemNumber)
for dp in range(len(patch.readings[rec].Vdp)):
    node1_x = patch.readings[rec].Vdp[dp].Rx1East
    node1_y = patch.readings[rec].Vdp[dp].Rx1North
    node2_x = patch.readings[rec].Vdp[dp].Rx2East
    node2_y = patch.readings[rec].Vdp[dp].Rx2North
    # plt.plot(node1_x, node1_y, 'o');
    # plt.plot(node2_x, node2_y, 'o');
    plt.plot([node1_x, node2_x], [node1_y, node2_y], '-o');

plt.title(str(patch.readings[rec].MemNumber))
plt.show()
