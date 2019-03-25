import matplotlib.pyplot as plt
import numpy as np
import DCIPtools as DCIP
from multiprocessing import Pool

# ============================ User area ==============================
# ver 0.1
# set variables
fname = "C:\\Users\\johnk\\Projects\\Karst\\Lusha_North-QC.DAT"
# load the data file
patch = DCIP.loadDias(fname)

# get a node database for a reading
dipole_angles = patch.get3Dangles(bin_width=20,
                                  dipole_max=175,
                                  dipole_min=15)
theta = np.asarray(dipole_angles)
# print(theta.size)
# Compute pie slices
bin_width = 10
bins_start = np.arange(0, 360, bin_width)
bins_end = np.arange(20, 380, bin_width)
print(bins_start)
print(bins_end)
theta = np.asarray(patch.readings[0].angles)
radii = np.zeros(bins_end.size)
for i in range(bins_start.size):
    accept = np.logical_and(theta > bins_start[i], theta < bins_end[i])
    accepted = theta[accept]
    radii[i] = accepted.size

ax = plt.subplot(111, projection='polar')
bars = ax.bar(bins_start * np.pi / 180, radii, width=np.deg2rad(bin_width), bottom=0.0)
ax.set_thetagrids(np.arange(0, 360, 20), labels=np.arange(0, 360, 20))
ax.set_theta_zero_location('N')
ax.set_theta_direction(-1)

plt.show()
