import sys
from PyQt5 import QtCore
import matplotlib.pyplot as plt
import numpy as np
import JDataObject as Jdata
from matplotlib.colors import LinearSegmentedColormap
from PyQt5.QtWidgets import QApplication

# ==============================================================
# Global variable to declare
vmin = 10.                                  # vmin_rho
vmax = 4000.                                 # vmax_rho
vsmin = 1.                                  # vmin_rho
vsmax = 20.                                 # vmax_rho
fig, ((ax, ax2), (ax3, ax4)) = plt.subplots(2, 2)
update = True
line_direction = 'north'
# Head to bottom of script to identify file path
# =============================================================


def file_changed(path):
    """
    method for updating plot when file changes

    """
    print('File Changed: %s' % path)
    #  ============================================================
    # define the file required for import
    fileName = path
    ax.clear()
    ax2.clear()
    ax3.clear()
    ax4.clear()
    injection_location_x = []
    injection_location_y = []
    xp_rho = []
    yp_rho = []
    xp_mx = []
    yp_mx = []
    val = []
    val2 = []
    id1_rho = []
    id2_rho = []
    id1_mx = []
    id2_mx = []
    xp_rho_rev = []
    yp_rho_rev = []
    xp_mx_rev = []
    yp_mx_rev = []
    val_rev = []
    val2_rev = []
    id1_rho_rev = []
    id2_rho_rev = []
    id1_mx_rev = []
    id2_mx_rev = []
    injection = []
    # =================================================================
    # Code Start
    patch = Jdata.loadDias(fileName)               # loads data
    print(len(patch.readings))
    # calculated mid-pt data points
    for src in range(len(patch.readings)):
        injection_location_x.append(patch.readings[src].Idp.Tx1)
        injection_location_y.append(25.0)
        injection.append(patch.readings[src].MemNumber)
        for rx in range(len(patch.readings[src].Vdp)):
            if line_direction == 'east':
                patch.readings[src].Vdp[rx].direction = False

            if patch.readings[src].Vdp[rx].flagRho == "Accept":
                try:
                    if (patch.readings[src].Idp.Tx1 <
                            patch.readings[src].Vdp[rx].Rx1):
                        xp_rho.append(
                            patch.readings[src].Vdp[rx].getXplotpoint(
                                patch.readings[src].Idp))
                        yp_rho.append(
                            patch.readings[src].Vdp[rx].getZplotpoint(
                                patch.readings[src].Idp))
                        val.append(
                            np.abs(patch.readings[src].Vdp[rx].Rho))
                        id1_rho.append(
                            patch.readings[src].Vdp[rx].Rx1File[:2])
                        id2_rho.append(
                            patch.readings[src].Vdp[rx].Rx2File[:2])
                    elif (patch.readings[src].Idp.Tx1 >
                            patch.readings[src].Vdp[rx].Rx1):
                        xp_rho_rev.append(
                            patch.readings[src].Vdp[rx].getXplotpoint(
                                patch.readings[src].Idp))
                        yp_rho_rev.append(
                            patch.readings[src].Vdp[rx].getZplotpoint(
                                patch.readings[src].Idp))
                        val_rev.append(
                            np.abs(patch.readings[src].Vdp[rx].Rho))
                        id1_rho_rev.append(
                            patch.readings[src].Vdp[rx].Rx1File[:2])
                        id2_rho_rev.append(
                            patch.readings[src].Vdp[rx].Rx2File[:2])
                except:
                    print("something wrong with file loading")
                    pass
            if patch.readings[src].Vdp[rx].flagMx == "Accept":
                try:
                    if (patch.readings[src].Idp.Tx1 <
                            patch.readings[src].Vdp[rx].Rx1):
                        xp_mx.append(
                            patch.readings[src].Vdp[rx].getXplotpoint(
                                patch.readings[src].Idp))
                        yp_mx.append(
                            patch.readings[src].Vdp[rx].getZplotpoint(
                                patch.readings[src].Idp))
                        val2.append(
                            np.abs(patch.readings[src].Vdp[rx].Mx))
                        id1_mx.append(
                            patch.readings[src].Vdp[rx].Rx1File[:2])
                        id2_mx.append(
                            patch.readings[src].Vdp[rx].Rx2File[:2])
                    elif (patch.readings[src].Idp.Tx1 >
                            patch.readings[src].Vdp[rx].Rx1):
                        xp_mx_rev.append(
                            patch.readings[src].Vdp[rx].getXplotpoint(
                                patch.readings[src].Idp))
                        yp_mx_rev.append(
                            patch.readings[src].Vdp[rx].getZplotpoint(
                                patch.readings[src].Idp))
                        val2_rev.append(
                            np.abs(patch.readings[src].Vdp[rx].Mx))
                        id1_mx_rev.append(
                            patch.readings[src].Vdp[rx].Rx1File[:2])
                        id2_mx_rev.append(
                            patch.readings[src].Vdp[rx].Rx2File[:2])
                except:
                    print("something wrong with file loading")
                    pass

    # convert to numpy
    midx = np.asarray(xp_rho)
    midz = np.asarray(yp_rho)
    midx_rev = np.asarray(xp_rho_rev)
    midz_rev = np.asarray(yp_rho_rev)
    injection_location_y = np.asarray(injection_location_y)
    injection_location_x = np.asarray(injection_location_x)
    # check which data to plot
    name = 'custom_div_cmap'
    pcolorOpts = {}
    custom_map = LinearSegmentedColormap.from_list(name=name,
                                                   colors=['aqua',
                                                           [0, 0.85, 1, 1],
                                                           [0.1, 1, 0.1, 1],
                                                           'yellow',
                                                           [1, 0.7, 0, 1],
                                                           [1, 0.2, 0.2, 1],
                                                           [0.95, 0.9, 1, 1]],
                                                   N=200)
    # ==========================================================
    # Plot forward direction Resistivity
    ax.axes.set_title("Apparent Resistivity (ohm-m)", y=1.14)
    try:
        ph = ax.scatter(xp_rho, yp_rho, c=val, cmap=custom_map,
                        vmin=vmin, vmax=vmax, clim=(vmin, vmax))
        #  plotting dipole info text
        for i, txt in enumerate(val):
            idizzle = id1_rho[i] + "-" + id2_rho[i]
            ax.annotate(idizzle, (midx[i], midz[i]), size=6)
        # Initiate color bar if needed
        if update:
            cbar = plt.colorbar(ph, ax=ax,
                                format="%.0f", fraction=0.04,
                                orientation="vertical")

    except ValueError:
        print("no data to plot in forward Rho")
    # =============================================================
    # Plot reverse direction Resistivity
    try:
        ph_reverse = ax3.scatter(xp_rho_rev, yp_rho_rev, c=val_rev,
                                 cmap=custom_map,
                                 vmin=vmin, vmax=vmax, clim=(vmin, vmax))
        #  plotting dipole info text
        for i, txt in enumerate(val_rev):
            idizzle2 = id1_rho_rev[i] + "-" + id2_rho_rev[i]
            ax3.annotate(idizzle2, (midx_rev[i], midz_rev[i]), size=6)
        # Initiate color bar if needed
        if update:
            cbar3 = plt.colorbar(ph_reverse, ax=ax3,
                                 format="%.0f", fraction=0.04,
                                 orientation="vertical")

    except ValueError:
        print("no data to plot in reverse Rho")
    #  ===========================================================
    #  plot the injection text
    for i in range(len(patch.readings)):
        ax.annotate(int(injection[i]), (injection_location_x[i],
                    injection_location_y[i]), size=6)
        ax2.annotate(int(injection[i]), (injection_location_x[i],
                     injection_location_y[i]), size=6)
        ax3.annotate(int(injection[i]), (injection_location_x[i],
                     injection_location_y[i]), size=6)
        ax4.annotate(int(injection[i]), (injection_location_x[i],
                     injection_location_y[i]), size=6)

    #  ===========================================================
    #  Plot forward direction Chargeability
    name = 'custom_div_cmap1'
    ax2.axes.set_title("Apparent Chargeability (mV/V)", y=1.11)
    try:
        ph2 = ax2.scatter(xp_mx, yp_mx, c=val2, cmap=custom_map,
                          vmin=vsmin, vmax=vsmax, clim=(vmin, vmax))
        #  plot dipole info text
        for i, txt in enumerate(val2):
            idizzle = id1_mx[i] + "-" + id2_mx[i]
            ax2.annotate(idizzle, (midx[i], midz[i]), size=6)
        # Initiate colorbar if needed
        if update:
            cbar2 = plt.colorbar(ph2, ax=ax2,
                                 format="%.0f", fraction=0.04,
                                 orientation="vertical")

    except ValueError:
            print("no data to plot in forward Mx")

    #  =========================================================
    #  Plot reverse direction Chargeability
    try:
        ph4 = ax4.scatter(xp_mx_rev, yp_mx_rev, c=val2_rev, cmap=custom_map,
                          vmin=vsmin, vmax=vsmax, clim=(vmin, vmax))
        #  plot dipole info text
        for i, txt in enumerate(val2_rev):
            idizzle2 = id1_mx_rev[i] + "-" + id2_mx_rev[i]
            ax4.annotate(idizzle2, (midx_rev[i], midz_rev[i]), size=6)
        # Initiate colorbar if needed
        if update:
            cbar4 = plt.colorbar(ph4, ax=ax4,
                                 format="%.0f", fraction=0.04,
                                 orientation="vertical")

    except ValueError:
            print("no data to plot in reverse Mx")

    # ============================================================
    #  adjust plotting axis for new data as it comes in
    try:
        ax.axes.get_yaxis().set_visible(False)
        ax.axes.set_ylim(np.min(midz) - 50, np.max(midz) + 100)
        ax.axes.set_xlim(np.min(midx) - 250, np.max(midx) + 100)
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')
        ax.tick_params(labelsize=8)
        ax2.axes.get_yaxis().set_visible(False)
        ax2.axes.set_ylim(np.min(midz) - 50, np.max(midz) + 100)
        ax2.axes.set_xlim(np.min(midx) - 250, np.max(midx) + 100)
        ax2.xaxis.tick_top()
        ax2.xaxis.set_label_position('top')
    except ValueError:
        pass
    try:
        ax2.tick_params(labelsize=8)
        ax3.axes.get_yaxis().set_visible(False)
        ax3.axes.set_ylim(np.min(midz_rev) - 50, np.max(midz_rev) + 100)
        ax3.axes.set_xlim(np.min(midx_rev) - 250, np.max(midx_rev) + 100)
        ax3.xaxis.tick_top()
        ax3.xaxis.set_label_position('top')
        ax3.tick_params(labelsize=8)
        ax4.axes.get_yaxis().set_visible(False)
        ax4.axes.set_ylim(np.min(midz_rev) - 50, np.max(midz_rev) + 100)
        ax4.axes.set_xlim(np.min(midx_rev) - 250, np.max(midx_rev) + 100)
        ax4.xaxis.tick_top()
        ax4.xaxis.set_label_position('top')
        ax4.tick_params(labelsize=8)
    except ValueError:
        pass
    plt.draw()
# ================================================================


# ==================================================================
# Running the Application
app = QApplication(sys.argv)

paths = ["C:\\Users\\johnk\\devProjects\\Python\\sectionplots\\L4250.dat"]
# "C:\Users\johnk\devProjects\Python\sectionplots\L0L1A.dat"
try:
    file_changed(paths[0])
    update = False
except:
    print("no file ready for plotting")

fs_watcher = QtCore.QFileSystemWatcher(paths)
# fs_watcher.directoryChanged.connect(directory_changed)
fs_watcher.fileChanged.connect(file_changed)
plt.show()

sys.exit(app.exec_())
