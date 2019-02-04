from simpegEMIP.StretchedExponential import SEInvProblem, SESurvey
from SimPEG import *
import matplotlib.pyplot as plt
import numpy as np
import DCIPtools as DCIP


def fit_with_se(time, dobs, eta0=0.01, tau0=0.1, c0=0.5):
    siginf = 1.
    wires = Maps.Wires(('eta', 1), ('tau', 1), ('c', 1))
    taumap = Maps.ExpMap(nP=1)*wires.tau
    etamap = Maps.ExpMap(nP=1)*wires.eta
    cmap = Maps.ExpMap(nP=1)*wires.c    
    survey = SESurvey()
    survey.dobs = dobs
    m1D = Mesh.TensorMesh([np.ones(3)])
    prob = SEInvProblem(m1D, etaMap = etamap, tauMap = taumap, cMap=cmap)
    prob.time = time
    prob.pair(survey)
    eta0 = dobs[1]
    m0 = np.log(np.r_[eta0, tau0, c0])
    perc = 0.05
    dmisfitpeta = DataMisfit.l2_DataMisfit(survey)
    dmisfitpeta.W = 1/(abs(survey.dobs)*perc)
    reg = Regularization.Simple(m1D)
    opt = Optimization.ProjectedGNCG(maxIter = 20)
    invProb = InvProblem.BaseInvProblem(dmisfitpeta, reg, opt)
    target = Directives.TargetMisfit()
    betaSch = Directives.BetaSchedule(coolingFactor=1, coolingRate=1)
    opt.upper = np.log(np.r_[1., 10., 1.])
    opt.lower = np.log(np.r_[1e-5, 0.001, 1e-2])
    
    invProb.beta=0.
    inv = Inversion.BaseInversion(invProb, directiveList=[target])
    prob.counter = opt.counter = Utils.Counter()
    opt.LSshorten = 0.5
    opt.remember('xc')
    opt.tolX = 1e-20
    opt.tolF = 1e-20
    opt.tolG = 1e-20
    opt.eps = 1e-20
    mopt = inv.run(m0)   
    return np.exp(mopt), survey.dobs, invProb.dpred


# ============================ User area ==============================
# set variables
fname = "/Users/juan/Documents/testData/L14_200-1200m_QC.DAT"
outname = "/Users/juan/Documents/testData/L14_200-1200m_cole-cole.DAT"
# load the data file
patch = DCIP.loadDias(fname)
# set the start and end time of the data that is fed into inversion
start_time = 400
end_time = 1760
# make time vector
times_gates = patch.window_center
# geometric center equition for time window
times = np.sqrt(patch.window_start * patch.window_end)
# get indicies
tinds = times > start_time
error = []
c = []
tau = []
eta = []
index_list = []
# run the inversion scheme
for rdg in range(len(patch.readings)):
    for dp in range(len(patch.readings[rdg].Vdp)):
        # if patch.readings[rdg].Vdp[dp].flagMx == "Accept":
        # try:
        Vp = patch.readings[rdg].Vdp[dp].Vp
        Vs = patch.readings[rdg].Vdp[dp].Vs
        times_gates = patch.window_center
        mopt, dobs, dpred = fit_with_se(times[tinds] * 1e-3, Vs[tinds] / Vp)
        mx_pred = np.sum(dpred * patch.window_width[tinds]) / np.sum(patch.window_width[tinds]) * 1e3
        mx_obs = np.sum(dobs * patch.window_width[tinds]) / np.sum(patch.window_width[tinds]) * 1e3
        patch.readings[rdg].Vdp[dp].Mx = mx_obs
        err = np.sqrt(np.sum((dobs - dpred)**2)) / dobs.size * 1e3
        np.abs(err / mx_obs) * 100
        print("eta: {0}, c: {1}, tau: {2} error eta: {3} : Mx: {4}".format(mx_pred, mopt[1], mopt[2], err, mx_obs))
        
        if np.isnan(err):
            error.append(-99.9)
            c.append(-99.9)
            tau.append(-99.9)
            eta.append(-99.9)
        else:
            error.append(err)
            c.append(mopt[1])
            tau.append(mopt[2])
            eta.append(mx_pred)
        index = np.array([rdg, dp])
        index_list.append(index)
        # except:
        #     pass

# create organised numpy arrays
c_ = np.vstack(c)
eta_ = np.vstack(eta)
tau_ = np.vstack(tau)
active_inds = np.vstack(index_list)
error = np.vstack(error)
plt.plot(error, '*')
plt.show()
# Figure out the times
start_inds = (patch.window_start > start_time)
stop_inds = (patch.window_end < end_time)
strt_time = patch.window_start[start_inds]
stop_time = patch.window_end[stop_inds]
# write the data to file
patch.writeColeColeDat(eta_, tau_, c_, error, active_inds, outname,
                       strt_time[0], stop_time[stop_time.size - 1])
# get app. res. data and plot it
Rho = patch.getApparentResistivity(reject="Mx")
fig, axs = plt.subplots(4, 1)
properties = [Rho, eta_,
              tau_, c_]
titles = ["$\\rho_{a}$", "$\\eta_{a}$", "$\\tau_{a}$", "$c_{a}$"]
colors = ['#1f77b4', 'seagreen', 'crimson', 'gold']
for i, ax in enumerate(axs):
    out = ax.hist(np.log10(properties[i]), bins=50, color=colors[i])
    ax.set_title(titles[i])
    ax.set_xticklabels([("%.1f") % (10**tick) for tick in ax.get_xticks()])
plt.tight_layout()
plt.show()
