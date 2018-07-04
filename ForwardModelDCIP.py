"""
3D DC forward model & inversion of pole-dipolearray
======================================
"""

from SimPEG import (
    Mesh, Maps, Utils,
    DataMisfit, Regularization, Optimization,
    InvProblem, Directives, Inversion
)
from SimPEG.EM.Static import DC, Utils as DCUtils
import numpy as np
import matplotlib.pyplot as plt
from pymatsolver import Pardiso as Solver
from time import clock
np.random.seed(12345)

# =============================
# methods


def getRxData():
    xyz = open("C:/Users/johnk/Projects/Seabridge/IdealizedStations_Rx.csv")
    x = []
    y = []
    z = []
    for line in xyz:
        x_, y_, z_ = line.split(',')
        x1_ = float(x_)
        y1_ = float(y_)
        dist = np.sqrt((x1_ - 374725.)**2 +
                       (y1_ - 6276216.)**2)
        if dist > 160.:
            x.append(float(x_))
            y.append(float(y_))
            z.append(float(z_))
    xyz.close()
    x1 = np.asarray(x)
    y1 = np.asarray(y)
    z1 = np.asarray(z)
    rx_electrodes = np.c_[x1, y1, z1]
    return rx_electrodes


def getTxData():
    xyz = open("C:/Users/johnk/Projects/Seabridge/IdealizedStations_Tx.csv")
    x = []
    y = []
    z = []
    for line in xyz:
        x_, y_, z_ = line.split(',')
        x1_ = float(x_)
        y1_ = float(y_)
        dist = np.sqrt((x1_ - 374725.)**2 +
                       (y1_ - 6276216.)**2)
        if dist > 160.:
            x.append(float(x_))
            y.append(float(y_))
            z.append(float(z_))
    xyz.close()
    x.append(374739.0)
    y.append(6276261.0)
    z.append(1855.0)
    x1 = np.asarray(x)
    y1 = np.asarray(y)
    z1 = np.asarray(z)
    tx_electrodes = np.c_[x1, y1, z1]
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
                distE = np.abs(node1[0] - tx[idk, 0])
                if distE < 650:
                    if (min_dipole_size) < dist < (max_dipole_size):
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


# ============================================================================
# start script
fileName1 = "C:/Users/johnk/Projects/Seabridge/fmdata.con"        # output mod
fileName2 = "C:/Users/johnk/Projects/Seabridge/forwardmodel.msh"  # input mesh
mesh = Mesh.TensorMesh._readUBC_3DMesh(fileName2)    # Read in/create the mesh

# =============================================================================
# create sphere for ice representation
x0 = (np.max(mesh.gridCC[:, 0]) +
      np.min(mesh.gridCC[:, 0])) / 2. + 50        # x0 center point of sphere
y0 = (np.max(mesh.gridCC[:, 1]) +
      np.min(mesh.gridCC[:, 1])) / 2. - 50        # y0 center point of sphere
z0 = 2350                                          # x0 center point of sphere
# (np.max(mesh.gridCC[:, 2]) + np.min(mesh.gridCC[:, 2])) / 2.
r0 = 500                                           # radius of sphere
print(x0, y0, z0)
csph = (np.sqrt((mesh.gridCC[:, 0] - x0)**2. +
                (mesh.gridCC[:, 1] - y0)**2. +
                (mesh.gridCC[:, 2] - z0)**2.)) < r0  # indicies of sphere
# sphere done =================================================================

print("Starting forward modeling")
start = clock()
# Define model Background
ln_sigback = np.log(1. / 33000.)                     # log conductivity
rx = getRxData()                                     # rx locations
tx = getTxData()                                     # tx locations
survey = generateSurvey(rx, tx, 45, 135)             # create survey object
survey.getABMN_locations()                           # get locations
uniq = Utils.uniqueRows(np.vstack((survey.a_locations,
                                   survey.b_locations,
                                   survey.m_locations,
                                   survey.n_locations)))  # row the locations
electrode_locations = uniq[0]                            # assign
actinds = Utils.surface2ind_topo(mesh,
                                 electrode_locations,
                                 method='cubic')      # active indicies
survey.drapeTopo(mesh, actinds)                       # drape topo
# ============================================================================
# Create model
sigma = np.ones(mesh.nC) * 1. / 33000.                # conductivity

# create dipping structure parameters
theta = 45. * np.pi / 180.                            # dipping angle
x0_d = 374400.
x1_d = 374950.
y0_d = 6276000.
y0_1d = 900. * np.sin(theta) + y0_d
y1_d = 6276100.
y1_1d = 900. * np.sin(theta) + y1_d
z0_d = 1710.
z1_d = z0_d - (700. * np.cos(theta))
m_ = (z0_d - z1_d) / (y0_1d - y0_d)                   # slope of dip

# loop through mesh and assign dipping structure conductivity
for idx in range(mesh.nC):
    if z1_d <= mesh.gridCC[idx, 2] <= z0_d:
        if (x0_d <= mesh.gridCC[idx, 0] <= x1_d):
            zslope1 = z0_d - m_ * (mesh.gridCC[idx, 1] - y0_d)
            zslope2 = z0_d - m_ * (mesh.gridCC[idx, 1] - y1_d)
            yslope1 = y0_d + (1. / m_) * (mesh.gridCC[idx, 2] - z0_d)
            yslope2 = y1_d + (1. / m_) * (mesh.gridCC[idx, 2] - z0_d)
            if yslope1 <= mesh.gridCC[idx, 1] <= yslope2:
                sigma[idx] = 1. / 300.

# sigma[csph] = ((1. / 5000.) *
#                np.ones_like(sigma[csph]))             # set sphere values
sigma[~actinds] = 1. / 1e8                            # flag air values
rho = 1. / sigma
stop = clock()
print(stop)
# plot results
ncy = mesh.nCy
ncz = mesh.nCz
ncx = mesh.nCx
mtrue = rho
clim = [0, 60000.]
fig, ax = plt.subplots(2, 2, figsize=(12, 6))
ax = Utils.mkvc(ax)
dat = mesh.plotSlice(((mtrue)), ax=ax[0], normal='Z', clim=clim,
                     ind=int(ncz / 2 - 4))
ax[0].plot(rx[:, 0], rx[:, 1], 'or')
ax[0].plot(tx[:, 0], tx[:, 1], 'dk')
mesh.plotSlice(((mtrue)), ax=ax[1], normal='Y', clim=clim,
               ind=int(ncy / 2))
mesh.plotSlice(((mtrue)), ax=ax[2], normal='X', clim=clim,
               ind=int(ncx / 2 + 4))
cbar_ax = fig.add_axes([0.82, 0.15, 0.05, 0.7])
cb = plt.colorbar(dat[0], ax=cbar_ax)
fig.subplots_adjust(right=0.85)
cb.set_label('rho')
cbar_ax.axis('off')
plt.show()

# End model creation =========================================================
# ============================================================================
actmap = Maps.InjectActiveCells(
    mesh, indActive=actinds,
    valInactive=np.log(1e8))                           # active map
mapping = Maps.ExpMap(mesh) * actmap                   # create the mapping
mtrue = np.log(sigma[actinds])                         # create the true model
problem = DC.Problem3D_CC(mesh, sigmaMap=mapping)      # assign problem
problem.pair(survey)                                   # pair the survey + prob
problem.Solver = Solver                                # assign solver
survey.dpred(mtrue)                                    # view predicted data
survey.makeSyntheticData(mtrue, std=0.05, force=True)  # make synthetic data
stop = clock()
print("timing of making synthetic data:", stop)
print("starting inversion of forward data")
# Tikhonov Inversion
####################
start = clock()
# Initial Model
m0 = np.median(ln_sigback) * np.ones(mapping.nP)
# Data Misfit
dmis = DataMisfit.l2_DataMisfit(survey)
uncert = abs(survey.dobs) * 0.05 * (10**(-3.2))
dmis.W = 1. / uncert
# Regularization
regT = Regularization.Tikhonov(mesh, indActive=actinds, alpha_s=1e-4,
                               alpha_x=1., alpha_y=1., alpha_z=1.)

# Optimization Scheme
opt = Optimization.InexactGaussNewton(maxIter=5)

# Form the problem
opt.remember('xc')
invProb = InvProblem.BaseInvProblem(dmis, regT, opt)

# Directives for Inversions
beta = Directives.BetaEstimate_ByEig(beta0_ratio=1e0)
Target = Directives.TargetMisfit()
betaSched = Directives.BetaSchedule(coolingFactor=5., coolingRate=2)

inv = Inversion.BaseInversion(invProb, directiveList=[beta, Target,
                                                      betaSched])
# Run Inversion
minv = inv.run(m0)
stop = clock()
print("timing of inversion:", stop)
out_model = np.ones(sigma.size)
out_model[~actinds] = 1. / 1e8
out_model[actinds] = np.exp(minv)
# write data to file
out_file = open(fileName1, "w")
for i in range(out_model.size):
    out_file.write("%0.5e\n" % out_model[i])

# plot results
ncy = mesh.nCy
ncz = mesh.nCz
ncx = mesh.nCx
mtrue = out_model
# print(mtrue.min(), mtrue.max())
clim = [0, 50000.]
fig, ax = plt.subplots(2, 2, figsize=(12, 6))
ax = Utils.mkvc(ax)
dat = mesh.plotSlice(((mtrue)), ax=ax[0], normal='Z', clim=clim,
                     ind=int(40))
ax[0].plot(rx[:, 0], rx[:, 1], 'or')
mesh.plotSlice(((mtrue)), ax=ax[1], normal='Y', clim=clim,
               ind=int(ncy / 2))
mesh.plotSlice(((mtrue)), ax=ax[2], normal='X', clim=clim,
               ind=int(ncx / 2))
cbar_ax = fig.add_axes([0.82, 0.15, 0.05, 0.7])
cb = plt.colorbar(dat[0], ax=cbar_ax)
fig.subplots_adjust(right=0.85)
cb.set_label('rho')
cbar_ax.axis('off')
plt.show()
