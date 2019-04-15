import matplotlib.pyplot as plt
import numpy as np
import DCIPtools as DCIP
import pylab as plt
from scipy.interpolate import griddata

################################################################
# define the file required for import
fileName = "F:/EMX/Deliverable/Pseudosections/Boggy_L8400_Final_50m2_A.DAT"
unitType = "appResistivity"
# unitType = "appChareability"
pd_config = True
xp = []
yp = []
val = []
val2 = []
maxN = 38
Nlevel = ["N = 1", "N = 2", "N = 3", "N = 4", "N = 5", "N = 6", "N = 7", "N = 8", "N = 9", "N = 10", "N = 11","N = 12", "N = 13"]
#, "N = 14", "N = 11","N = 12", "N = 13", "N = 14", "N = 15", "N = 16", "N = 17", "N = 18", "N = 19", "N = 20", "N = 21", "N = 22", "N = 23", "N = 24", "N = 25", "N = 26", "N = 27", "N = 28", "N = 29", "N = 30", "N = 31", "N = 32", "N = 33", "N = 34", "N = 35", "N = 36", "N = 37", "N = 38"]
#, "N = 14", "N = 15", "N = 16", "N = 17", "N = 18", "N = 19", "N = 20", "N = 21", "N = 22", "N = 23", "N = 24", "N = 25", "N = 26", "N = 27", "N = 28", "N = 29", "N = 30", "N = 31", "N = 32", "N = 33", "N = 34", "N = 35", "N = 36", "N = 37", "N = 38", "N = 39", "N = 40", "N = 41", "N = 42", "N = 43", "N = 44", "N = 45", "N = 46", "N = 47", "N = 48", "N = 49", "N = 50", "N = 51", "N = 52"
# z_n = [-50, -75, -100, -125, -150, -175, -200]
z_n = -25.0 * np.arange(1, maxN + 1) - 50
vmin_rho, vmax_rho = 5, 800
vmin_mx, vmax_mx = 0, 20

# =================================================================
# Code Start
patch = DCIP.loadDias(fileName)               # loads data
# calculated mid-pt data points
for src in range(len(patch.readings)):
    for rx in range(len(patch.readings[src].Vdp)):
        if patch.readings[src].Vdp[rx].flagMx == "Accept":
            xp.append(
                patch.readings[src].Vdp[rx].getXplotpoint(patch.readings[src].Idp, dipole_dipole=pd_config))
            yp.append(
                patch.readings[src].Vdp[rx].getZplotpoint(patch.readings[src].Idp, dipole_dipole=pd_config))
            val.append(
                np.abs(patch.readings[src].Vdp[rx].Rho))
            val2.append(
                np.abs(patch.readings[src].Vdp[rx].Mx))
# convert to numpy
midx = np.asarray(xp)
midz = np.asarray(yp)
rho = np.asarray(val)
mx = np.asarray(val2)
xNLevel = np.min(midx) - 200
x_n = np.zeros(len(z_n))
x_n = x_n + xNLevel
# Grid points
grid_x, grid_z = np.mgrid[np.min(midx):np.max(midx),
                          np.min(midz):np.max(midz)]
# create an axes to plot on
ax = plt.subplot(2, 1, 1, aspect='equal')
# create an axes to plot on
ax2 = plt.subplot(2, 1, 2, aspect='equal')
# check which data to plot
# if unitType == "appResistivity":
vmin = vmin_rho
vmax = vmax_rho
ax.axes.set_title("Apparent Resistivity (ohm-m)", y=1.14)
grid_rho = griddata(np.c_[midx, midz], rho.T, (grid_x, grid_z),
                    method='cubic')
grid_rho = np.ma.masked_where(np.isnan(grid_rho), grid_rho)
name = 'custom_div_cmap'
pcolorOpts = {}
from matplotlib.colors import LinearSegmentedColormap
custom_map = LinearSegmentedColormap.from_list(name=name, colors=['aqua', [0, 0.85, 1, 1], [0.1, 1, 0.1, 1], 'yellow', [1, 0.7, 0, 1], [1, 0.2, 0.2, 1], [0.95, 0.9, 1, 1]], N=200)
CS = ax.contour(grid_x[:, 0], grid_z[0, :], grid_rho.T, 15, linewidths=0.5, colors='k')
ph = ax.pcolormesh(grid_x[:, 0], grid_z[0, :], grid_rho.T, cmap=custom_map, clim=(vmin, vmax), vmin=vmin, vmax=vmax, **pcolorOpts)
# , vmin=vmin, vmax=vmax, {})
cbar = plt.colorbar(ph, ax=ax, format="%.0f", fraction=0.04, orientation="vertical")
cbar.set_label("App.Res.", size=8)
for i, txt in enumerate(val):
    ax.annotate(int(rho[i]), (midx[i], midz[i]), size=7)
# else:
vmin = vmin_mx
vmax = vmax_mx
ax2.axes.set_title("Apparent Chargeability (mV/V)", y=1.11)
grid_mx = griddata(np.c_[midx, midz], mx.T, (grid_x, grid_z),
                   method='cubic')
grid_mx = np.ma.masked_where(np.isnan(grid_mx), grid_mx)
name = 'custom_div_cmap1'
pcolorOpts = {}
# from matplotlib.colors import LinearSegmentedColormap
# custom_map = LinearSegmentedColormap.from_list(name=name, colors=['aqua', [0, 0.85, 1, 1], [0.1, 1, 0.1, 1], 'yellow', [1, 0.7, 0, 1], [1, 0.2, 0.2, 1], [0.95, 0.9, 1, 1]], N=200)
CS2 = ax2.contour(grid_x[:, 0], grid_z[0, :], grid_mx.T, 15, linewidths=0.5, colors='k')
ph2 = ax2.pcolormesh(grid_x[:, 0], grid_z[0, :], grid_mx.T, cmap=custom_map, clim=(vmin, vmax), vmin=vmin, vmax=vmax, **pcolorOpts)
# , vmin=vmin, vmax=vmax, {})
cbar2 = plt.colorbar(ph2, ax=ax2, format="%.0f", fraction=0.04, orientation="vertical")
cbar2.set_label("App.Mx.", size=8)
for i, txt in enumerate(val):
    ax2.annotate(int(mx[i]), (midx[i], midz[i]), size=7)

for i, txt in enumerate(Nlevel):
    ax.annotate(Nlevel[i], (x_n[i], z_n[i]), size=8)
    ax2.annotate(Nlevel[i], (x_n[i], z_n[i]), size=8)

ax.axes.get_yaxis().set_visible(False)
ax.axes.set_ylim(np.min(midz) - 50, np.max(midz) + 50)
ax.axes.set_xlim(np.min(midx) - 250, np.max(midx) + 100)
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')
ax.tick_params(labelsize=8)
ax2.axes.get_yaxis().set_visible(False)
ax2.axes.set_ylim(np.min(midz) - 50, np.max(midz) + 50)
ax2.axes.set_xlim(np.min(midx) - 250, np.max(midx) + 100)
ax2.xaxis.tick_top()
ax2.xaxis.set_label_position('top')
ax2.tick_params(labelsize=8)
plt.show()
