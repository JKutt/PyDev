from SimPEG import DC, IP
from SimPEG import Maps, Utils
import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
from pylab import hist
try:
    from pymatsolver import Pardiso as Solver
except ImportError:
    from SimPEG import SolverLU as Solver


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
                if (min_dipole_size) < dist < (max_dipole_size):
                    if node1[0] != tx[idk, 0]:
                        if node1[0] < tx[idk, 0] < node2[0]:
                            print("inbetween")
                        else:
                            rx1.append(node1)
                            rx2.append(node2)
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
    background = 15000.
    # Obtain ABMN locations
    # Generate DC survey object
    rx1 = np.arange(100, 450, 50)
    rx2 = np.arange(750, 1050, 50)
    rx = np.r_[rx1, rx2]
    # rx = np.arange(100, 1000, 50)
    rx_y = np.zeros(rx.size)
    rx = np.c_[rx, rx_y]
    tx1 = np.arange(50, 500, 50)
    tx2 = np.arange(700, 1100, 50)
    tx3 = np.ones(1) * 600
    tx = np.r_[tx1, tx2]
    # tx = np.arange(50, 1050, 50)
    tx_y = np.zeros(tx.size)
    tx = np.c_[tx, tx_y]
    survey_dc = generateSurvey(rx, tx, 45, 55)
    # survey_dc = DC.Utils.gen_DCIPsurvey(endl, survey_type=survey_type, dim=2,
    #                                     a=50, b=50, n=18)
    # survey_dc.survey_type = survey_type
    survey_dc.getABMN_locations()
    survey_dc = IO.from_ambn_locations_to_survey(
        survey_dc.a_locations, survey_dc.b_locations,
        survey_dc.m_locations, survey_dc.n_locations,
        survey_type, data_dc_type='volt', data_ip_type='volt'
    )

    # Obtain 2D TensorMesh
    mesh, actind = IO.set_mesh()
    topo, mesh1D = DC.Utils.genTopography(mesh, -50, 0, its=100)
    fill = topo[0:topo.size / 3]  # messing with topo
    topo[topo.size / 3:(topo.size / 3 + topo.size / 3)] = fill  # messing with topo
    actind = Utils.surface2ind_topo(mesh, np.c_[mesh1D.vectorCCx, topo])
    survey_dc.drapeTopo(mesh, actind, option="top")

    # Build conductivity and chargeability model ==========
    blk_inds_charg = Utils.ModelBuilder.getIndicesSphere(
        np.r_[600., -80.], 12.5, mesh.gridCC
    )
    blk_inds_charg2 = Utils.ModelBuilder.getIndicesSphere(
        np.r_[580., -90.], 12.5, mesh.gridCC
    )
    blk_inds_charg3 = Utils.ModelBuilder.getIndicesSphere(
        np.r_[570., -100.], 12.5, mesh.gridCC
    )
    blk_inds_charg4 = Utils.ModelBuilder.getIndicesSphere(
        np.r_[560., -110.], 12.5, mesh.gridCC
    )
    blk_inds_charg5 = Utils.ModelBuilder.getIndicesSphere(
        np.r_[550., -120.], 12.5, mesh.gridCC
    )
    blk_inds_charg6 = Utils.ModelBuilder.getIndicesSphere(
        np.r_[540., -130.], 12.5, mesh.gridCC
    )
    blk_inds_charg7 = Utils.ModelBuilder.getIndicesSphere(
        np.r_[530., -140.], 12.5, mesh.gridCC
    )
    blk_inds_charg8 = Utils.ModelBuilder.getIndicesSphere(
        np.r_[520., -150.], 12.5, mesh.gridCC
    )
    blk_inds_charg9 = Utils.ModelBuilder.getIndicesSphere(
        np.r_[510., -160.], 12.5, mesh.gridCC
    )
    blk_inds_charg10 = Utils.ModelBuilder.getIndicesSphere(
        np.r_[490., -170.], 12.5, mesh.gridCC
    )
    blk_inds_charg11 = Utils.ModelBuilder.getIndicesSphere(
        np.r_[470., -180.], 12.5, mesh.gridCC
    )
    # conductivity model
    sigma = np.ones(mesh.nC) * 1./background
    sigma[~actind] = 1./1e8
    sigma[blk_inds_charg] = 1. / 100.
    sigma[blk_inds_charg2] = 1. / 100.
    sigma[blk_inds_charg3] = 1. / 100.
    sigma[blk_inds_charg4] = 1. / 100.
    sigma[blk_inds_charg5] = 1. / 100.
    sigma[blk_inds_charg6] = 1. / 100.
    sigma[blk_inds_charg7] = 1. / 100.
    sigma[blk_inds_charg8] = 1. / 100.
    sigma[blk_inds_charg9] = 1. / 100.
    sigma[blk_inds_charg10] = 1. / 100.
    sigma[blk_inds_charg11] = 1. / 100.
    rho = 1./sigma    # resistivity
    # chargeability
    charg = np.ones(mesh.nC) * 0.018
    charg[blk_inds_charg] = 0.035
    charg[blk_inds_charg2] = 0.035
    charg[blk_inds_charg3] = 0.035
    charg[blk_inds_charg4] = 0.035
    charg[blk_inds_charg5] = 0.035
    charg[blk_inds_charg6] = 0.035
    charg[blk_inds_charg7] = 0.035
    charg[blk_inds_charg8] = 0.035
    charg[blk_inds_charg9] = 0.035
    charg[blk_inds_charg10] = 0.035
    charg[blk_inds_charg11] = 0.035
    # End target creation =========================

    # Show the true conductivity model
    if plotIt:
        fig, axs = plt.subplots(2,1, figsize=(12, 6))
        temp_rho = rho.copy()
        temp_rho[~actind] = np.nan
        temp_charg = charg.copy()
        temp_charg[~actind] = np.nan

        out1 = mesh.plotImage(
            temp_rho, grid=True, ax=axs[0], gridOpts={'alpha': 0.2},
            clim=(10, 1000),
            pcolorOpts={"cmap": "jet", "norm": colors.LogNorm()}
        )
        out2 = mesh.plotImage(
            temp_charg, grid=True, ax=axs[1], gridOpts={'alpha': 0.2},
            clim=(0, 0.1),
            pcolorOpts={"cmap": "jet"}
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

    # Use Exponential Map: m = log(rho)
    actmap = Maps.InjectActiveCells(
        mesh, indActive=actind, valInactive=np.log(1e8)
    )
    mapping = Maps.ExpMap(mesh) * actmap

    # Generate mtrue_dc for resistivity
    mtrue_dc = np.log(rho[actind])

    # Generate 2D DC problem
    prb = DC.Problem2D_N(
        mesh, rhoMap=mapping, storeJ=True,
        Solver=Solver
    )
    # Pair problem with survey
    try:
        prb.pair(survey_dc)
    except:
        survey_dc.unpair()
        prb.pair(survey_dc)

    # Make synthetic DC data with 5% Gaussian noise
    dtrue_dc = survey_dc.makeSyntheticData(mtrue_dc, std=0.06, force=True)
    IO.data_dc = dtrue_dc

    # Generate mtrue_ip for chargability
    mtrue_ip = charg[actind]
    # Generate 2.5D DC problem
    # "N" means potential is defined at nodes
    # prb_ip = IP.Problem2D_N(
    #     mesh, etaMap=actmap, storeJ=True, rho=rho,
    #     Solver=Solver
    # )
    survey_ip = IP.from_dc_to_ip_survey(survey_dc, dim="2D")
    # prb_ip.pair(survey_ip)
    # dtrue_ip = survey_ip.makeSyntheticData(mtrue_ip, std=0.06, force=True)
    # IO.data_ip = dtrue_ip

    # Show apparent resisitivty pseudo-section
    if plotIt:
        IO.plotPseudoSection(
            data_type='apparent_resistivity', scale='linear',
            cmap='jet'
        )
        plt.show()

    # Show apparent resisitivty histogram
    if plotIt:
        fig = plt.figure(figsize=(10, 4))
        ax1 = plt.subplot(121)
        out = hist(np.log10(abs(IO.voltages)), bins=20)
        ax1.set_xlabel("log10 DC voltage (V)")
        ax2 = plt.subplot(122)
        out = hist(IO.apparent_resistivity, bins=20)
        ax2.set_xlabel("Apparent Resistivity ($\Omega$m)")
        plt.tight_layout()
        plt.show()

    # Set initial model
    m0_dc = np.ones(actmap.nP) * np.log(background)
    # Set uncertainty
    # floor
    eps_dc = 10**(-4.2)
    # percentage
    std_dc = 0.06

    mopt_dc, pred_dc = DC.run_inversion(
        m0_dc, survey_dc, actind, mesh, std_dc, eps_dc,
        use_sensitivity_weight=False
        )
    # plt.hist(dtrue_dc)
    # plt.show()
    # DC.Utils.writeUBC_DCobs("fmdc.obs", survey_dc, 2, 'SIMPLE', 'pole-dipole')
    # outmodel = 1./ rho
    # outmodel[actind] = mopt_dc
    # outmodel[~actind] = 1e-8
    # mesh.writeUBC("fm.msh")
    # mesh.writeModelUBC("cond-fm.con", outmodel)
    # Convert obtained inversion model to resistivity
    # rho = M(m), where M(.) is a mapping
    rho_est = mapping*mopt_dc
    # rho_est[~actind] = np.nan
    rho_true = rho.copy()
    rho_true[~actind] = np.nan

    # show recovered conductivity
    if plotIt:
        vmin, vmax = rho.min(), rho.max()
        fig, ax = plt.subplots(2, 1, figsize=(20, 6))
        out1 = mesh.plotImage(
                rho_true, clim=(10, 15000),
                pcolorOpts={"cmap": "jet"},
                ax=ax[0]
        )
        out2 = mesh.plotImage(
            rho_est, clim=(10, 15000),
            pcolorOpts={"cmap": "jet"},
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

    # DC.Utils.writeUBC_DCobs("fmip.obs", survey_ip, 2, 'GENERAL', 'pole-dipole')
    prb_ip = IP.Problem2D_N(
        mesh, etaMap=actmap, storeJ=False, rho=rho,
        Solver=Solver
    )
    prb_ip.pair(survey_ip)
    dtrue_ip = survey_ip.makeSyntheticData(mtrue_ip, std=0.06, force=True)
    IO.data_ip = dtrue_ip

    # Show apparent chargeability pseudo-section
    if plotIt:
        IO.plotPseudoSection(
            data_type='apparent_chargeability', scale='linear',
            cmap='jet'
        )
        plt.show()

    # Show apparent Mx histogram
    if plotIt:
        fig = plt.figure(figsize=(10, 4))
        ax1 = plt.subplot(121)
        out = hist(np.log10(abs(IO.voltages_ip)), bins=20)
        ax1.set_xlabel("log10 IP voltage (V)")
        ax2 = plt.subplot(122)
        out = hist(IO.apparent_chargeability, bins=20)
        ax2.set_xlabel("Apparent Chargeability (V/V)")
        plt.tight_layout()
        plt.show()
    # Set initial model based upon histogram
    m0_ip = np.ones(actmap.nP) * 0.018
    # Set uncertainty
    # floor
    eps_ip = 10**(-4)
    # percentage
    std_ip = 0.06
    # Clean sensitivity function formed with true resistivity
    prb_ip._Jmatrix = None
    # Input obtained resistivity to form sensitivity
    prb_ip.rho = mapping * mopt_dc
    mopt_ip, _ = IP.run_inversion(
        m0_ip, survey_ip, actind, mesh, std_ip, eps_ip,
        upper=np.Inf, lower=0.,
        beta0_ratio=1e2,
        use_sensitivity_weight=False
    )

    # Convert obtained inversion model to chargeability
    # charg = M(m), where M(.) is a mapping for cells below topography

    charg_est = actmap*mopt_ip
    charg_est[~actind] = np.nan
    charg_true = charg.copy()
    charg_true[~actind] = np.nan

    # show recovered chargeability
    if plotIt:
        fig, ax = plt.subplots(2, 1, figsize=(20, 6))
        out1 = mesh.plotImage(
                charg_true, clim=(0, 0.04),
                pcolorOpts={"cmap": "jet"},
                ax=ax[0]
        )
        out2 = mesh.plotImage(
            charg_est, clim=(0, 0.04),
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
            cb.set_label("Mx (V/V)")
            ax[i].set_xlabel("Northing (m)")
            ax[i].set_ylabel("Elevation (m)")
            ax[i].set_aspect('equal')
        plt.tight_layout()
        plt.show()

if __name__ == '__main__':
    survey_type = 'pole-dipole'
    run(survey_type=survey_type, plotIt=True)
