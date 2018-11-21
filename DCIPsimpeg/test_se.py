from simpegEMIP.StretchedExponential import SEInvProblem, SESurvey
from SimPEG import *
import matplotlib.pyplot as plt
import numpy as np
import DCIPtools as DCIP
# from SimPEG import DC


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


fname = "/Users/juan/Documents/testData/Seabridge_Final.DAT"
outname = "/Users/juan/Documents/testData/Seabridge-Final-cole-cole.DAT"
format_dias = '%s,%s,%s,%0.3e,%0.3e,%0.3e,%0.3e\n' 
patch = DCIP.loadDias(fname)
out_file = open(outname, "w+")
out_file.write('%s\n' % ("record,id1,id2,c,tau,eta,err"))
rdg = 0
dp = 4

Vp = patch.readings[rdg].Vdp[dp].Vp
Vs = patch.readings[rdg].Vdp[dp].Vs
times_gates = patch.window_center
times = np.sqrt(patch.window_start * patch.window_end) * 1e-3
print(Vp)
print(Vs / Vp)
# plt.plot(times, Vs, 'o-')
mopt, dobs, dpred = fit_with_se(times, Vs / Vp * -1)
mx_pred = np.sum(dpred * patch.window_width) / np.sum(patch.window_width) * 1e3
err = np.sqrt(np.sum((dobs - dpred)**2)) / dobs.size * 1e3
print("eta: {0}, c: {1}, tau: {2} error eta: {3}".format(mx_pred, mopt[1], mopt[2], err))

# out_file.write(format_dias % (patch.readings[rdg].Vdp[dp].reading,
#                               patch.readings[rdg].Vdp[dp].Rx1File[:2],
#                               patch.readings[rdg].Vdp[dp].Rx2File[:2],
#                               mopt[1], mopt[2], mx_pred, err))
plt.plot(times, dobs, '*')
plt.plot(times, dpred, 'r')
plt.show()
# error = []
# c = []
# tau = []
# eta_ = []
# for rdg in range(1):
#     for dp in range(len(patch.readings[rdg].Vdp)):
#         try:
#             Vp = patch.readings[rdg].Vdp[dp].Vp
#             Vs = patch.readings[rdg].Vdp[dp].Vs
#             times_gates = patch.window_center
#             times = np.sqrt(patch.window_start * patch.window_end) * 1e-3
#             # print(Vp)
#             # print(Vs[6:] / Vp)
#             # plt.plot(times[6:], Vs[6:], 'o-')
#             mopt, dobs, dpred = fit_with_se(times, Vs / Vp)
#             mx_pred = np.sum(dpred * patch.window_width) / np.sum(patch.window_width) * 1e3
#             err = np.sqrt(np.sum((dobs - dpred)**2)) / dobs.size * 1e3
#             print("eta: {0}, c: {1}, tau: {2} error eta: {3}".format(mx_pred, mopt[1], mopt[2], err))
#             c.append(mopt[1])
#             tau.append(mopt[2])
#             eta_.append(mx_pred)
#             error.append(err)
#             out_file.write(format_dias % (patch.readings[rdg].Vdp[dp].reading,
#                                           patch.readings[rdg].Vdp[dp].Rx1File[:2],
#                                           patch.readings[rdg].Vdp[dp].Rx2File[:2],
#                                           mopt[1], mopt[2], mx_pred, err))
#         except:
#             pass
# error = np.asarray(error)
# c = np.asarray(c)
# eta_ = np.asarray(eta_)
# tau = np.asarray(tau)
# # plt.hist(error)
# # plt.show()

# plt.scatter(c, tau, c=eta_,
#             cmap='viridis')
# plt.colorbar()  # show color scale
# plt.show()
