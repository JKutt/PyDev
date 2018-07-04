
from SimPEG import DC, IP
from SimPEG import Maps, Utils
from SimPEG import Mesh
import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
from time import clock
from pylab import hist
try:
    from pymatsolver import Pardiso as Solver
except ImportError:
    from SimPEG import SolverLU as Solver


def getRxData():
    xyz = open("/Users/juan/Documents/testData/IdealizedStations_Rx.csv")
    x = []
    y = []
    z = []
    for line in xyz:
        x_, y_, z_ = line.split(',')
        x1_ = float(x_)
        y1_ = float(y_)
        dist = np.sqrt((x1_ - 374775.)**2 +
                       (y1_ - 6276016.)**2)

        if 250 > dist > 10:
            x.append(float(x_))
            y.append(float(y_))
            z.append(float(z_))
            # z.append(2176.)
    xyz.close()
    x1 = np.asarray(x)
    y1 = np.asarray(y)
    z1 = np.asarray(z)
    rx_electrodes = np.c_[x1, y1, z1]
    return rx_electrodes


def getTxData():
    xyz = open("/Users/juan/Documents/testData/IdealizedStations_Tx.csv")
    x = []
    y = []
    z = []
    for line in xyz:
        x_, y_, z_ = line.split(',')
        x1_ = float(x_)
        y1_ = float(y_)
        dist = np.sqrt((x1_ - 374775.)**2 +
                       (y1_ - 6276016.)**2)
        if 250 > dist > 10:
            x.append(float(x_))
            y.append(float(y_))
            z.append(float(z_))
            # z.append(2176.)
        # x.append(float(x_))
        # y.append(float(y_))
        # z.append(float(z_))
    xyz.close()
    # x.append(374739.0)
    # y.append(6276261.0)
    # z.append(1855.0)
    x1 = np.asarray(x)
    y1 = np.asarray(y)
    z1 = np.asarray(z)
    tx_electrodes = np.c_[x1, y1, z1]
    tx_electrodes = tx_electrodes[::2]
    return tx_electrodes


def generateSurvey(rx, tx, min_dipole_size, max_dipole_size):
    """
     Generates a survey to through into a forward model
     INPUT:
          rx_dx = array of Rx x spacings
          rx_dy = array of Rx y spacings
          Tx_dx = array of Tx x spacings
          Tx_dy = array of Tx y spacings
    """
    SrcList = []
    rx_length = rx.shape[0]
    for idk in range(tx.shape[0]):
        rx1 = []
        rx2 = []
        for idx in range(rx_length):
            node1 = rx[idx, :]
            for idj in range(idx, rx_length):
                node2 = rx[idj, :]
                dist = np.sqrt(np.sum((node1 - node2)**2))
                distr1 = np.sqrt(np.sum((node1 - tx[idk, :])**2))
                distr2 = np.sqrt(np.sum((node2 - tx[idk, :])**2))
                distE = np.abs(node1[0] - tx[idk, 0])
                if distE < 80:
                    if (min_dipole_size) < dist < (max_dipole_size):
                        if distr1 < distr2:
                            rx1.append(node1)
                            rx2.append(node2)
                        else:
                            rx1.append(node2)
                            rx2.append(node1)
                    # print(dist)
        rx1 = np.asarray(rx1)
        rx2 = np.asarray(rx2)
        rxClass = DC.Rx.Dipole(rx1, rx2)
        srcClass = DC.Src.Pole([rxClass], tx[idk, :])
        SrcList.append(srcClass)

    survey = DC.Survey(SrcList)

    return survey


def run(plotIt=True, survey_type="pole-dipole"):
    np.random.seed(1)
    # Initiate I/O class for DC
    IO = DC.IO()
    # Obtain ABMN locations

    fileName1 = "/Users/juan/Documents/testData/fmdata-daf.con"        # output mod
    fileName2 = "/Users/juan/Documents/testData/forwardmodel.msh"  # input mesh
    mesh = Mesh.TensorMesh._readUBC_3DMesh(fileName2)     # Read in/create mesh

    print("Starting forward modeling")
    start = clock()
    # Define model Background
    rx = getRxData()                                     # rx locations
    tx = getTxData()                                     # tx locations
    survey_dc = generateSurvey(rx, tx, 45, 55)         # create survey object
    survey_dc.survey_type = "pole-dipole"
    survey_dc.getABMN_locations()                           # get locations
    # survey_dc = IO.from_ambn_locations_to_survey(
    #     survey_dc.a_locations, survey_dc.b_locations,
    #     survey_dc.m_locations, survey_dc.n_locations,
    #     survey_type, data_dc_type='volt', data_ip_type='volt'
    # )
    uniq = Utils.uniqueRows(np.vstack((survey_dc.a_locations,
                                       survey_dc.b_locations,
                                       survey_dc.m_locations,
                                       survey_dc.n_locations)))
    electrode_locations = uniq[0]                            # assign
    actinds = Utils.surface2ind_topo(mesh,
                                     electrode_locations,
                                     method='cubic')      # active indicies
    # survey_dc.drapeTopo(mesh, actinds)
    IO.a_locations = survey_dc.a_locations.copy()
    IO.b_locations = survey_dc.b_locations.copy()
    IO.m_locations = survey_dc.m_locations.copy()
    IO.n_locations = survey_dc.n_locations.copy()          # drape topo
    IO.data_dc_type = 'volt'
    IO.G = IO.geometric_factor(survey_dc)
# =============================================================================
    # create sphere for ice representation
    x0 = (np.max(mesh.gridCC[:, 0]) +
          np.min(mesh.gridCC[:, 0])) / 2. + 50      # x0 center point of sphere
    y0 = (np.max(mesh.gridCC[:, 1]) +
          np.min(mesh.gridCC[:, 1])) / 2. - 50      # y0 center point of sphere
    z0 = 2350                                       # x0 center point of sphere
    # (np.max(mesh.gridCC[:, 2]) + np.min(mesh.gridCC[:, 2])) / 2.
    r0 = 500                                           # radius of sphere
    print(x0, y0, z0)
    csph = (np.sqrt((mesh.gridCC[:, 0] - x0)**2. +
                    (mesh.gridCC[:, 1] - y0)**2. +
                    (mesh.gridCC[:, 2] - z0)**2.)) < r0  # indicies of sphere
# sphere done =================================================================
# ============================================================================
# Create model
    mx = np.ones(mesh.nC) * 0.018                         # chargeability
    sigma = np.ones(mesh.nC) * 1. / 15000.

# create dipping structure parameters
    theta = 45. * np.pi / 180.                            # dipping angle
    x0_d = 374700.
    x1_d = 375000.
    y0_d = 6275850.
    y0_1d = 900. * np.sin(theta) + y0_d
    y1_d = 6275900.
    y1_1d = 900. * np.sin(theta) + y1_d
    z0_d = 1900.
    z1_d = z0_d - (900. * np.cos(theta))
    m_ = (z0_d - z1_d) / (y0_1d - y0_d)                   # slope of dip

    # loop through mesh and assign dipping structure conductivity
    for idx in range(mesh.nC):
        if z1_d <= mesh.gridCC[idx, 2] <= z0_d:
            if (x0_d <= mesh.gridCC[idx, 0] <= x1_d):
                yslope1 = y0_d + (1. / m_) * (mesh.gridCC[idx, 2] - z0_d)
                yslope2 = y1_d + (1. / m_) * (mesh.gridCC[idx, 2] - z0_d)
                if yslope1 <= mesh.gridCC[idx, 1] <= yslope2:
                    mx[idx] = 0.035
                    sigma[idx] = 1. / 300.

    # mx[csph] = ((0.025) *
    #             np.ones_like(mx[csph]))             # set sphere values
    mx[~actinds] = 1. / 1e8                             # flag air values
    # sigma[csph] = ((5000.) *
    #                np.ones_like(sigma[csph]))             # set sphere values
    sigma[~actinds] = 1. / 1e8                              # flag air values
    rho = 1. / sigma
    stop = clock()
    print(stop)
    # plot results
    # Show the true conductivity model
    if plotIt:
        ncy = mesh.nCy
        ncz = mesh.nCz
        ncx = mesh.nCx
        print(mesh.nC)
        clim = [0, 0.04]
        fig, ax = plt.subplots(2, 2, figsize=(12, 6))
        ax = Utils.mkvc(ax)
        dat = mesh.plotSlice(((mx)), ax=ax[0], normal='Z', clim=clim,
                             ind=int(ncz / 2 - 14), pcolorOpts={"cmap": "jet"})
        ax[0].plot(rx[:, 0], rx[:, 1], 'or')
        ax[0].plot(tx[:, 0], tx[:, 1], 'dk')
        mesh.plotSlice(((mx)), ax=ax[1], normal='Y', clim=clim,
                       ind=int(ncy / 2 + 2), pcolorOpts={"cmap": "jet"})
        mesh.plotSlice(((mx)), ax=ax[2], normal='X', clim=clim,
                       ind=int(ncx / 2 + 4), pcolorOpts={"cmap": "jet"})
        mesh.plotSlice(((mx)), ax=ax[3], normal='X', clim=clim,
                       ind=int(ncx / 2 + 8), pcolorOpts={"cmap": "jet"})
        cbar_ax = fig.add_axes([0.82, 0.15, 0.05, 0.7])
        cb = plt.colorbar(dat[0], ax=cbar_ax)
        fig.subplots_adjust(right=0.85)
        cb.set_label('V/V')
        cbar_ax.axis('off')
        plt.show()
        # print(mtrue.min(), mtrue.max())
        # clim = [0, 20000]
        # fig, ax = plt.subplots(2, 2, figsize=(12, 6))
        # ax = Utils.mkvc(ax)
        # dat = mesh.plotSlice(((mtrue)), ax=ax[0], normal='Z', clim=clim,
        #                      ind=int(ncz / 2 - 4), pcolorOpts={"cmap": "jet"})
        # ax[0].plot(rx[:, 0], rx[:, 1], 'or')
        # mesh.plotSlice(((mtrue)), ax=ax[1], normal='Y', clim=clim,
        #                ind=int(ncy / 2), pcolorOpts={"cmap": "jet"})
        # mesh.plotSlice(((mtrue)), ax=ax[2], normal='X', clim=clim,
        #                ind=int(ncx / 2 + 4), pcolorOpts={"cmap": "jet"})
        # mesh.plotSlice(((mtrue)), ax=ax[3], normal='X', clim=clim,
        #                ind=int(ncx / 2 + 8), pcolorOpts={"cmap": "jet"})
        # cbar_ax = fig.add_axes([0.82, 0.15, 0.05, 0.7])
        # cb = plt.colorbar(dat[0], ax=cbar_ax)
        # fig.subplots_adjust(right=0.85)
        # cb.set_label('rho')
        # cbar_ax.axis('off')
        # plt.show()

    # print(error.size, survey_ip.obs.size)
    # error = np.asarray(survey_ip.dobs) * (np.random.randint(-1, 2, survey_ip.dobs) / 20.)
    # survey_ip.dobs = survey_ip.dobs + error
    # Use Exponential Map: m = log(rho)
    actmap = Maps.InjectActiveCells(
        mesh, indActive=actinds, valInactive=np.log(1e8)
    )
    mapping = Maps.ExpMap(mesh) * actmap

    # Generate mtrue_dc for resistivity
    mtrue_dc = np.log(rho[actinds])

    # Generate 3D DC problem
    # "CC" means potential is defined at center
    prb = DC.Problem3D_CC(
        mesh, rhoMap=mapping, storeJ=False,
        Solver=Solver
    )
    # Pair problem with survey
    # try:
    prb.pair(survey_dc)
    # except:
    #     survey_dc.unpair()
    #     prb.pair(survey_dc)

    survey_dc.dpred(mtrue_dc)
    # Make synthetic DC data with 5% Gaussian noise
    dtrue_dc = survey_dc.makeSyntheticData(mtrue_dc,
                                           std=0.05,
                                           force=True)
    IO.data_dc = dtrue_dc
    # Generate mtrue_ip for chargability
    mtrue_ip = mx[actinds]
    # Generate 3D DC problem
    # "CC" means potential is defined at center
    prb_ip = IP.Problem3D_CC(
        mesh, etaMap=actmap, storeJ=False, rho=rho,
        Solver=Solver
    )
    survey_ip = IP.from_dc_to_ip_survey(survey_dc, dim="3D")
    survey_ip.survey_type = "pole-dipole"
    survey_ip.getABMN_locations()
    prb_ip.pair(survey_ip)
    survey_ip.dpred(mtrue_ip)
    dtrue_ip = survey_ip.makeSyntheticData(mtrue_ip, std=0.05)
    survey_ip.std = np.ones_like(dtrue_ip) * 0.05
    survey_ip.eps = np.ones_like(dtrue_ip) * 10**(-4)
    IO.data_ip = dtrue_ip
    DC.Utils.writeUBC_DCobs("fmip.obs", survey_ip, 3, 'GENERAL', 'pole-dipole')
    # # Show apparent resisitivty histogram
    # if plotIt:
        # plt.plot(dtrue_dc, 'o')
        # fig = plt.figure(figsize=(10, 4))
        # ax1 = plt.subplot(121)
        # out = hist(np.log10(abs(IO.voltages)), bins=20)
        # ax1.set_xlabel("log10 DC voltage (V)")
        # ax2 = plt.subplot(122)
        # out = hist(IO.apparent_resistivity, bins=20)
        # ax2.set_xlabel("Apparent Resistivity ($\Omega$m)")
        # plt.tight_layout()
        # plt.show()


if __name__ == '__main__':
    survey_type = 'pole-dipole'
    run(survey_type=survey_type, plotIt=True)
