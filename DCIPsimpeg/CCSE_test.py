from simpegEMIP.StretchedExponential import SEInvProblem, SESurvey
from SimPEG import *
import matplotlib.pyplot as plt
import numpy as np
import DCIPtools as DCIP


def fit_with_se(time, dobs, eta0=0.01, tau0=0.1, c0=0.5):
    siginf = 1.
    wires = Maps.Wires(('eta', 1), ('tau', 1), ('c', 1))
    taumap = Maps.ExpMap(nP=1) * wires.tau
    etamap = Maps.ExpMap(nP=1) * wires.eta
    cmap = Maps.ExpMap(nP=1) * wires.c
    survey = SESurvey()
    survey.dobs = dobs
    m1D = Mesh.TensorMesh([np.ones(3)])
    prob = SEInvProblem(m1D, etaMap=etamap, tauMap=taumap, cMap=cmap)
    prob.time = time
    prob.pair(survey)
    eta0 = dobs[1]
    m0 = np.log(np.r_[eta0, tau0, c0])
    perc = 0.05
    dmisfitpeta = DataMisfit.l2_DataMisfit(survey)
    dmisfitpeta.W = 1 / (abs(survey.dobs) * perc)
    reg = Regularization.Simple(m1D)
    opt = Optimization.ProjectedGNCG(maxIter=20)
    invProb = InvProblem.BaseInvProblem(dmisfitpeta, reg, opt)
    target = Directives.TargetMisfit()
    betaSch = Directives.BetaSchedule(coolingFactor=1, coolingRate=1)
    opt.upper = np.log(np.r_[1., 10., 1.])
    opt.lower = np.log(np.r_[1e-5, 0.001, 1e-2])
    invProb.beta = 0.
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

if __name__ == '__main__':
    # ============================ User area ==============================
    # set variables
    fname = "E:/Projects/debug/Gunjan/L14_200-1200snappedPyCC400-1760.DAT"
    outname = "/Users/juan/Documents/testData/L14_200-1200m_cole-cole.DAT"
    # load the data file
    patch = DCIP.loadDias(fname)
    # set the start and end time of the data that is fed into inversion
    start_time = 400
    end_time = 1760
    rdg = 20
    dp = 87
    # make time vector
    times_gates = patch.window_center
    # geometric center equition for time window
    times = np.sqrt(patch.window_start * patch.window_end)
    # get indicies
    tinds = times > start_time
    Vp = patch.readings[rdg].Vdp[dp].Vp
    Vs = patch.readings[rdg].Vdp[dp].Vs
    times_gates = patch.window_center
    mopt, dobs, dpred = fit_with_se(times[tinds] * 1e-3, Vs[tinds] / Vp)
    mx_pred = np.sum(dpred * patch.window_width[tinds]) / np.sum(patch.window_width[tinds]) * 1e3
    mx_obs = np.sum(dobs * patch.window_width[tinds]) / np.sum(patch.window_width[tinds]) * 1e3
    err = np.sqrt(np.sum((dobs - dpred)**2)) / dobs.size * 1e3
    np.abs(err / mx_obs) * 100
    print("eta: {0}, c: {1}, tau: {2} error eta: {3} : Mx: {4}".format(mx_pred, mopt[1], mopt[2], err, mx_obs))
    print("previous error: {0}".format(patch.readings[rdg].Vdp[dp].Mx_err))
    print(patch.readings[rdg].Vdp[dp].reading, patch.readings[rdg].Vdp[dp].Rx1File)
    print(patch.readings[rdg].Vdp[dp].reading, patch.readings[rdg].Vdp[dp].Rx2File) 
    plt.plot(times[tinds], dobs, 'r')
    plt.plot(times[tinds], dpred, 'g')
    plt.show()
    