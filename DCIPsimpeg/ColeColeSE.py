import matplotlib.pyplot as plt
import numpy as np
import DCIPtools as DCIP
from simpegEMIP.StretchedExponential import SEMultiInvProblem, SEMultiSurvey
from SimPEG import (Mesh, Maps, Utils, DataMisfit, Regularization,
                    Optimization, Inversion, InvProblem, Directives)
from pymatsolver import PardisoSolver


def fitWithStretchedExponetial(time, survey_ip, sources, Rho, start_time=None):
    if start_time is None:
        start_time = 0
    # Choose all time channels
    tinds = time > start_time
    nLoc = survey_ip.dobs.T[tinds, :].shape[1]
    # Setup wire for different properties
    wires = Maps.Wires(('eta', nLoc), ('tau', nLoc), ('c', nLoc))
    taumap = Maps.ExpMap(nP=nLoc) * wires.tau
    etamap = Maps.ExpMap(nP=nLoc) * wires.eta
    cmap = Maps.ExpMap(nP=nLoc) * wires.c

    # This is almost dummmy mesh and xyz loc at the moment
    # But, there is potential use later
    m1D = Mesh.TensorMesh([np.ones(nLoc)])

    # # Set survey
    survey = SEMultiSurvey(time[tinds] * 1e-3, sources[:, :], n_pulse=2, T=4)
    survey.dobs = (survey_ip.dobs.T[tinds, :]).flatten(order='F')
    # Set problem
    prob = SEMultiInvProblem(m1D, etaMap=etamap, tauMap=taumap, cMap=cmap)
    prob.pair(survey)

    # Set initial model
    eta0, tau0, c0 = abs(survey_ip.dobs.T[0, :].T), 1. * np.ones(nLoc), 1. * np.ones(nLoc)

    m0 = np.r_[np.log(eta0), np.log(tau0), np.log(c0)]
    std = 0.02

    plt.plot(survey.dpred(m0))
    plt.plot(survey.dobs, '.')
    plt.title("Obs & pred")
    plt.xlabel("data point")
    plt.ylabel("value")
    plt.show()

    mreg = Mesh.TensorMesh([len(m0)])
    dmisfit = DataMisfit.l2_DataMisfit(survey)
    uncert = (abs(survey.dobs) * std + 0.01) * 1.1
    dmisfit.W = 1. / uncert
    reg = Regularization.Simple(mreg)
    opt = Optimization.ProjectedGNCG(maxIter=20)
    invProb = InvProblem.BaseInvProblem(dmisfit, reg, opt)
    # Create an inversion object
    target = Directives.TargetMisfit()
    invProb.beta = 0.
    inv = Inversion.BaseInversion(invProb, directiveList=[target])
    # reg.mref = 0.*m0
    prob.counter = opt.counter = Utils.Counter()
    opt.LSshorten = 0.5
    opt.remember('xc')
    opt.tolX = 1e-20
    opt.tolF = 1e-20
    opt.tolG = 1e-20
    opt.eps = 1e-20
    m_upper = np.r_[np.ones_like(tau0) * np.log(200),
                    np.ones_like(tau0) * np.log(10.),
                    np.ones_like(tau0) * np.log(1.)]
    m_lower = np.r_[np.ones_like(tau0) * np.log(1),
                    np.ones_like(tau0) * np.log(1e-2),
                    np.ones_like(tau0) * np.log(0.1)]

    opt.lower = m_lower
    opt.upper = m_upper
    mopt = inv.run(m0)

    eta = etamap * mopt
    tau = taumap * mopt
    c = cmap * mopt
    DPRED = invProb.dpred.reshape((nLoc, time[tinds].size))
    DOBS = survey.dobs.reshape((nLoc, time[tinds].size))
    UNCERT = uncert.reshape((nLoc, time[tinds].size))
    error = ((eta * 1e-3) / np.sqrt(np.sum((DOBS - DPRED)**2, axis=1) / time[tinds].size))

    from matplotlib import colors
    fig = plt.figure(figsize=(6, 5))
    active_inds = (Rho[:,0] > 10.) & (Rho[:,0] < 1e3) & (abs((DPRED - DOBS) / UNCERT).sum(axis=1) / time[tinds].size <1.)
    out = plt.scatter(tau[active_inds], c[active_inds],
                      c=eta[active_inds], norm=colors.LogNorm(),
                      cmap="magma", s=2)
    cb = plt.colorbar(out)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel("Tau")
    plt.ylabel("c")
    cb.set_label("Eta (mV/V)")
    plt.show()

    inds = (abs((DPRED - DOBS) / UNCERT).sum(axis=1) / time[tinds].size < 20.)
    # print(inds)
    inds = np.arange(survey.n_location)[inds]
    out = plt.loglog(time[tinds], DPRED[inds[::2], :].T, 'r')
    out = plt.loglog(time[tinds], DOBS[inds[::2], :].T, 'kx')
    print(np.sum((DPRED - DOBS)**2) / time.size)

    fig, axs = plt.subplots(4, 1)
    properties = [Rho[:, 0][active_inds], eta[active_inds],
                  tau[active_inds], c[active_inds]]
    titles = ["$\\rho_{a}$", "$\\eta_{a}$", "$\\tau_{a}$", "$c_{a}$"]
    colors = ['#1f77b4', 'seagreen', 'crimson', 'gold']
    for i, ax in enumerate(axs):
        out = ax.hist(np.log10(properties[i]), bins=50, color=colors[i])
        ax.set_title(titles[i])
        ax.set_xticklabels([("%.1f") % (10**tick) for tick in ax.get_xticks()])
    plt.tight_layout()
    plt.show()

    return eta, tau, c, error


# == Ver. 0.2
# ================================== User Area ==================
# File name
fname = "E:/Projects/debug/Charlotte/L4300_Final_All_NW.DAT"
# title of output file
outname = "E:/Projects/debug/Charlotte/L4300_Final_All_NW-CC.DAT"
# load the data file
patch = DCIP.loadDias(fname)
# create simpeg survey Object
survey_ip = patch.createDcSurvey("IP", ip_type="decay")
# create time vector
time = patch.window_center
# set the starting and ending times of the data for the inversion
start_time = 400.0                                                    #this one sets the start time for inversion of decay (e.g first window to start with)
end_time = 1760.0                                                     #and this is the last window to use.
# get all available injection locations from data file
injections = patch.getSources2(dipole=True)
# get all resistivities of accepted Mx data points
app_rho = patch.getApparentResistivity(reject="Mx")
# get indiiceis of all accepted Mx data points
inds_active = patch.getActiveIndicies(reject="Mx")
# run the inversion scheme
eta, tau, c, error = fitWithStretchedExponetial(time, survey_ip, injections,
                                                app_rho, start_time)
# write data to file
# patch.writeColeColeDat(eta, tau, c, error, inds_active, outname,
#                        start_time, end_time)
