from SimPEG import DC, IP
from SimPEG import Maps, Utils
from SimPEG import Mesh
import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
from time import clock
from pylab import hist
import DCIPtools as tools
try:
    from pymatsolver import Pardiso as Solver
except ImportError:
    from SimPEG import SolverLU as Solver


def run():
    np.random.seed(1)
    # Initiate I/O class for DC
    IO = DC.IO()
    IO_ip = DC.IO()
    filename = "/Users/juan/Documents/testData/L0_A.DAT"
    # filename1 = "/Users/juan/Documents/testData/KSA/ksa3D-DC.con"
    # filenameObs = "/Users/juan/Documents/testData/KSA/ksa3D-IP.obs"
    # fileName2 = "/Users/juan/Documents/testData/Mesh.msh"
    # mesh = Mesh.TensorMesh.readUBC(fileName2)
    # rho_est = mesh._readModelUBC_3D(filename1)
    patch = tools.loadDias(filename)
    survey_dc, tx = patch.createDcSurvey2D("DC")
    survey_ip, ipdat = patch.createDcSurvey2D("IP")
    survey_dc.dobs = tx
    survey_ip.dobs = ipdat
    IO.data_dc = tx
    IO.data_ip = ipdat
    IO_ip.data_ip = ipdat
    # plt.plot(survey_dc.dobs, '.')
    # for i in range(survey_dc.nSrc):
    #     for j in range(len(survey_dc.srcList[ii].rxList[0].locs)):
    #         rx = survey_dc.srcList[ii].rxList[0].locs
    #         rx1 = rx[:, 0]
    #         rx2 = rx[:, 1]
    #         # plt.plot(tx1[0], tx1[1], 'o')
    #         plt.plot([rx1[:, 0], rx2[:, 0]], [rx1[:, 1], rx2[:, 1]], '-o')

    # plt.show()
    # survey_dc = DC.Utils.readUBC_DC3Dobs(filenameObs)
    survey_type = "dipole-dipole"
    # survey_ip.data_dc_type = 'volt'
    survey_dc.getABMN_locations()
    survey_ip.getABMN_locations()
    # Obtain 2D TensorMesh
    uniq = Utils.uniqueRows(np.vstack((survey_dc.a_locations,
                                       survey_dc.b_locations,
                                       survey_dc.m_locations,
                                       survey_dc.n_locations)))
    electrode_locations = uniq[0]
    survey_dc = IO.from_ambn_locations_to_survey(
        survey_dc.a_locations, survey_dc.b_locations,
        survey_dc.m_locations, survey_dc.n_locations,
        survey_type, data_dc_type='volt',
        data_ip_type='apparent_chargeability'
    )
    survey_ip = IO_ip.from_ambn_locations_to_survey(
        survey_ip.a_locations, survey_ip.b_locations,
        survey_ip.m_locations, survey_ip.n_locations,
        survey_type, data_dc_type='volt',
        data_ip_type='apparent_chargeability'
    )
    survey_ip = IP.from_dc_to_ip_survey(survey_ip, dim="2.5D")
    survey_ip.dobs = ipdat
    mesh, actind = IO.set_mesh()
    # mesh.plotGrid()
    actinds = Utils.surface2ind_topo(mesh,
                                    electrode_locations,
                                    method='cubic')
    survey_dc.drapeTopo(mesh, actind)
    survey_ip.drapeTopo(mesh, actind)
    # IO.data_dc = survey_dc.dobs
    IO.plotPseudoSection(
        data_type='apparent_resistivity', scale='log',
        cmap='jet'
    )
    plt.show()

    # Show apparent chargeability pseudo-section

    IO_ip.plotPseudoSection(
        data_type='apparent_chargeability', scale='linear',
        cmap='jet'
    )
    plt.show()

    fig = plt.figure(figsize=(10, 4))
    ax1 = plt.subplot(121)
    out = hist(np.log10(abs(IO.voltages)), bins=20)
    ax1.set_xlabel("log10 DC voltage (V)")
    ax2 = plt.subplot(122)
    out = hist(IO.apparent_resistivity, bins=20)
    ax2.set_xlabel("Apparent Resistivity ($\Omega$m)")
    plt.tight_layout()
    plt.show()

    sigma = np.ones(mesh.nC) * 1. / 50.
    rho = 1. / sigma
    charg = np.ones(mesh.nC) * 0.002
    # Use Exponential Map: m = log(rho)
    actmap = Maps.InjectActiveCells(
        mesh, indActive=actind, valInactive=np.log(1e8)
    )
    mapping = Maps.ExpMap(mesh) * actmap

    # Generate mtrue_dc for resistivity
    mtrue_dc = np.log(rho[actind])

    fig, axs = plt.subplots(2, 1, figsize=(12, 6))
    temp_rho = rho.copy()
    temp_rho[~actind] = np.nan
    temp_charg = charg.copy()
    temp_charg[~actind] = np.nan

    out1 = mesh.plotImage(
        temp_rho, grid=True, ax=axs[0], gridOpts={'alpha': 0.5},
        clim=(10, 1000),
        pcolorOpts={"cmap": "viridis", "norm": colors.LogNorm()}
    )
    out2 = mesh.plotImage(
        temp_charg, grid=True, ax=axs[1], gridOpts={'alpha': 0.2},
        clim=(0, 0.1),
        pcolorOpts={"cmap": "magma"}
    )
    for i in range(2):
        axs[i].plot(
            survey_dc.electrode_locations[:, 0],
            survey_dc.electrode_locations[:, 1], 'kv'
        )
        axs[i].set_xlim(IO.grids[:, 0].min(), IO.grids[:, 0].max())
        axs[i].set_ylim(-IO.grids[:, 1].max(), IO.grids[:, 1].min())
        axs[i].set_aspect('equal')
    cb = plt.colorbar(out1[0], ax=axs[0])
    cb.set_label("Resistivity (ohm-m)")
    cb = plt.colorbar(out2[0], ax=axs[1])
    cb.set_label("Chargeability")

    plt.show()

    # Generate 2.5D DC problem
    # "N" means potential is defined at nodes
    prb = DC.Problem2D_N(
        mesh, rhoMap=mapping, storeJ=True,
        Solver=Solver
    )
    survey_dc.dobs = IO.data_dc
    prb.pair(survey_dc)

    # now for IP
    # Generate mtrue_ip for chargability
    mtrue_ip = charg[actind]
    # Generate 2.5D DC problem
    # "N" means potential is defined at nodes
    prb_ip = IP.Problem2D_N(
        mesh, etaMap=actmap, storeJ=True, rho=rho,
        Solver=Solver
    )
    survey_ip.dobs = IO_ip.data_ip
    prb_ip.pair(survey_ip)
    # Set initial model based upon histogram
    m0_dc = np.ones(actmap.nP) * np.log(50.)
    # Set uncertainty
    # floor
    eps_dc = 10**(-3.2)
    # percentage
    std_dc = 0.07

    mopt_dc, pred_dc = DC.run_inversion(
        m0_dc, survey_dc, actind, mesh, std_dc, eps_dc,
        beta0_ratio=1e0,
        use_sensitivity_weight=False
    )

    rho_est = mapping * mopt_dc
    rho_est[~actind] = np.nan
    rho_true = rho.copy()
    rho_true[~actind] = np.nan

    # show recovered conductivity
    plotIt = True
    if plotIt:
        vmin, vmax = rho.min(), rho.max()
        fig, ax = plt.subplots(2, 1, figsize=(30, 6))
        out1 = mesh.plotImage(
                rho_true, clim=(10, 400),
                pcolorOpts={"cmap": "jet", "norm": colors.LogNorm()},
                ax=ax[0]
        )
        out2 = mesh.plotImage(
            rho_est, clim=(10, 400),
            pcolorOpts={"cmap": "jet", "norm": colors.LogNorm()},
            ax=ax[1]
        )
        out = [out1, out2]
        for i in range(2):
            ax[i].plot(
                survey_dc.electrode_locations[:, 0],
                survey_dc.electrode_locations[:, 1], 'kv'
            )
            ax[i].set_xlim(IO.grids[:, 0].min(), IO.grids[:, 0].max())
            ax[i].set_ylim(-IO.grids[:, 1].max(), IO.grids[:, 1].min())
            cb = plt.colorbar(out[i][0], ax=ax[i])
            cb.set_label("Resistivity ($\Omega$m)")
            ax[i].set_xlabel("Northing (m)")
            ax[i].set_ylabel("Elevation (m)")
            ax[i].set_aspect('equal')
        plt.tight_layout()
        plt.show()

    # Set initial model based upon histogram
    m0_ip = np.ones(actmap.nP) * 1e-3
    # Set uncertainty
    # floor
    eps_ip = 10**(-4)
    # percentage
    std_ip = 0.07
    # Clean sensitivity function formed with true resistivity
    prb_ip._Jmatrix = None
    # Input obtained resistivity to form sensitivity
    prb_ip.rho = mapping * mopt_dc
    mopt_ip, _ = IP.run_inversion(
        m0_ip, survey_ip, actind, mesh, std_ip, eps_ip,
        maxIter=10,
        upper=np.Inf, lower=0.,
        beta0_ratio=1e0,
        use_sensitivity_weight=True,
    )

    # Convert obtained inversion model to chargeability
    # charg = M(m), where M(.) is a mapping for cells below topography

    charg_est = actmap * mopt_ip
    charg_est[~actind] = np.nan
    charg_true = charg.copy()
    charg_true[~actind] = np.nan

    # show recovered chargeability
    if plotIt:
        fig, ax = plt.subplots(2, 1, figsize=(20, 6))
        out1 = mesh.plotImage(
            charg_true, clim=(0, 0.06),
            pcolorOpts={"cmap": "jet"},
            ax=ax[0]
        )
        out2 = mesh.plotImage(
            charg_est, clim=(0, 0.06),
            pcolorOpts={"cmap": "jet"},
            ax=ax[1]
        )
        out = [out1, out2]
        for i in range(2):
            ax[i].plot(
                survey_dc.electrode_locations[:, 0],
                survey_dc.electrode_locations[:, 1], 'rv'
            )
            ax[i].set_xlim(IO.grids[:, 0].min(), IO.grids[:, 0].max())
            ax[i].set_ylim(-IO.grids[:, 1].max(), IO.grids[:, 1].min())
            cb = plt.colorbar(out[i][0], ax=ax[i])
            cb.set_label("Chargeability ($\Omega$m)")
            ax[i].set_xlabel("Northing (m)")
            ax[i].set_ylabel("Elevation (m)")
            ax[i].set_aspect('equal')
        plt.tight_layout()
        plt.show()


if __name__ == '__main__':
    survey_type = 'dipole-dipole'
    run()
