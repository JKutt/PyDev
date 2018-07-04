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

filename = "/Users/juan/Documents/testData/KSA/3D_QC_Preliminary.DAT"
filename1 = "/Users/juan/Documents/testData/KSA/ksa3D-DC.con"
filenameObs = "/Users/juan/Documents/testData/KSA/ksa3D-IP.obs"
fileName2 = "/Users/juan/Documents/testData/Mesh.msh"
mesh = Mesh.TensorMesh._readUBC_3DMesh(fileName2)
patch = tools.loadDias(filename)
print(len(patch.readings))
survey_dc, tx = patch.createDcSurvey("DC")
survey_dc.survey_type = "dipole-dipole"
survey_dc.data_dc_type = 'volt'
survey_dc.getABMN_locations()

# uniq = Utils.uniqueRows(np.vstack((survey_dc.a_locations,
#                                    survey_dc.b_locations,
#                                    survey_dc.m_locations,
#                                    survey_dc.n_locations)))
# electrode_locations = uniq[0]                            # assign
# actinds = Utils.surface2ind_topo(mesh,
#                                  electrode_locations,
#                                  method='cubic')      # active indicies
# survey_dc.drapeTopo(mesh, actinds)

survey_ip1, tx = patch.createDcSurvey("IP")
plt.plot(survey_ip1.dobs, '.')
plt.show()
# # Use Exponential Map: m = log(rho)
# actmap = Maps.InjectActiveCells(
#     mesh, indActive=actinds, valInactive=np.log(1e8)
# )
# mapping = Maps.ExpMap(mesh) * actmap

# # Generate 3D DC problem
# # "CC" means potential is defined at center
# prb = DC.Problem3D_CC(
#     mesh, rhoMap=mapping, storeJ=False,
#     Solver=Solver
# )
# # Pair problem with survey
# prb.pair(survey_dc)

# # Set initial model based upon histogram
# m0_dc = np.ones(actmap.nP) * np.log(406.)
# # Set uncertainty
# # floor
# eps_dc = 10**(-3.2)
# # percentage
# std_dc = 0.05
# mopt_dc, pred_dc = DC.run_inversion(
#     m0_dc, survey_dc, actinds, mesh, std_dc, eps_dc,
#     beta0_ratio=1e0,
#     use_sensitivity_weight=False
# )

# # Convert obtained inversion model to resistivity
# # rho = M(m), where M(.) is a mapping
# rho_est = mapping * mopt_dc
# # mesh.writeModelUBC(filename1, 1. / rho_est)
# mtrue = rho_est
# mtrue[~actinds] = np.nan
# clim = [0, 4000]
# ncy = mesh.nCy
# ncz = mesh.nCz
# ncx = mesh.nCx
# fig, ax = plt.subplots(2, 2, figsize=(12, 6))
# ax = Utils.mkvc(ax)
# dat = mesh.plotSlice(((mtrue)), ax=ax[0], normal='Z', clim=clim,
#                      ind=int(ncz / 2 + 6), pcolorOpts={"cmap": "jet"})
# ax[0].plot(tx[:, 0], tx[:, 1], 'or')
# mesh.plotSlice(((mtrue)), ax=ax[1], normal='Y', clim=clim,
#                ind=int(ncy / 2), pcolorOpts={"cmap": "jet"})
# # ax[0].plot(tx[:, 0], tx[:, 1], 'ok')
# mesh.plotSlice(((mtrue)), ax=ax[2], normal='X', clim=clim,
#                ind=int(ncx / 2 + 4), pcolorOpts={"cmap": "jet"})

# mesh.plotSlice(((mtrue)), ax=ax[3], normal='X', clim=clim,
#                ind=int(ncx / 2 + 10), pcolorOpts={"cmap": "jet"})
# cbar_ax = fig.add_axes([0.82, 0.15, 0.05, 0.7])
# cb = plt.colorbar(dat[0], ax=cbar_ax)
# fig.subplots_adjust(right=0.85)
# cb.set_label('rho')
# cbar_ax.axis('off')
# plt.show()

# # ip inversion =======================================
# prb_ip = IP.Problem3D_CC(
#     mesh, etaMap=actmap, storeJ=False, rho=rho_est,
#     Solver=Solver
# )
# survey_ip1, tx = patch.createDcSurvey("IP")
# survey_ip = IP.from_dc_to_ip_survey(survey_dc, dim="3D")
# # setting 
# survey_ip.survey_type = "dipole-dipole"
# survey_ip.dobs = survey_ip1.dobs
# survey_ip.std = survey_ip1.std
# survey_ip.eps = survey_ip1.eps
# survey_ip.getABMN_locations()
# survey_ip.drapeTopo(mesh, actinds)
# # survey_ip.dpred
# prb_ip.pair(survey_ip)

# # Set initial model based upon histogram
# m0_ip = np.ones(actmap.nP) * 1e-5
# # Set uncertainty
# # floor
# eps_ip = 10**(-3)
# # percentage
# std_ip = 0.05
# # Clean sensitivity function formed with true resistivity
# prb_ip._Jmatrix = None
# # Input obtained resistivity to form sensitivity
# prb_ip.rho = mapping * mopt_dc
# mopt_ip, _ = IP.run_inversion(
#     m0_ip, survey_ip, actinds, mesh, std_ip, eps_ip,
#     upper=np.Inf, lower=0.,
#     beta0_ratio=1e2,
#     use_sensitivity_weight=False
# )

# # Convert obtained inversion model to chargeability
# # charg = M(m), where M(.) is a mapping for cells below topography

# charg_est = actmap * mopt_ip
# mtrue = charg_est
# mtrue[~actinds] = np.nan
# clim = [0, 0.04]
# ncy = mesh.nCy
# ncz = mesh.nCz
# ncx = mesh.nCx
# fig, ax = plt.subplots(2, 2, figsize=(12, 6))
# ax = Utils.mkvc(ax)
# dat = mesh.plotSlice(((mtrue)), ax=ax[0], normal='Z', clim=clim,
#                      ind=int(ncz / 2 + 5), pcolorOpts={"cmap": "jet"})
# ax[0].plot(tx[:, 0], tx[:, 1], 'or')
# mesh.plotSlice(((mtrue)), ax=ax[1], normal='Y', clim=clim,
#                ind=int(ncy / 2), pcolorOpts={"cmap": "jet"})
# # ax[0].plot(tx[:, 0], tx[:, 1], 'ok')
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
