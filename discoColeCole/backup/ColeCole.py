# Discovery Cole-Cole script
# Imports ============================================
import JDataObject as jdata
import matplotlib.pyplot as plt
from time import clock
# load data ==========================================
fname = "C:\Projects\Dundee\proDATA\Deliverables\L3000E-PD.DAT"
outname = "C:\Projects\Dundee\proDATA\Deliverables\L3000E-PD-1_.DAT"
patch = jdata.loadDias(fname)
c_init = 0.3
tau_init = 0.000000025
radius = 1.0
# replace header
patch.headers[3] = "     RDG    DIPOLE     n     Status   Nominal_a         Tx1        Tx1East       Tx1North        Tx1Elev        Tx2East       Tx2North    Rx1FileName         Rx1        Rx1East       Rx1North        Rx1Elev   Rx2FileName          Rx2        Rx2East       Rx2North        Rx2Elev        Contact             SP             Vp         Vp_err             In         In_err            Rho        Rho_err         Rho_QC       Stack    TimeBase         MA     MA_err   MA_QC        C      Tau        M       M01       M02       M03       M04       M05       M06       M07       M08       M09       M10       M11       M12       M13       M14       M15       M16       M17       M18       M19       M20"
# end of load data ==================================

#  Cole-Cole fitting ================================
print("Starting cole-cole fitting")
start = clock()
for i in range(len(patch.readings)):
    for j in range(len(patch.readings[i].Vdp)):
        if patch.readings[i].Vdp[j].flagMx != "Reject":
            c, tau, M, err = (
                jdata.calcColeCole(patch.readings[i].Vdp[j].Vs /
                                   ((patch.readings[i].Vdp[j].Vp)),
                                   patch.window_width,
                                   c_init,
                                   tau_init,
                                   radius))
            patch.readings[i].Vdp[j].cole_c = c
            patch.readings[i].Vdp[j].cole_tau = tau
            patch.readings[i].Vdp[j].cole_m = M
            patch.readings[i].Vdp[j].Mx_err = err
            # print(c, tau, M, err / patch.readings[i].Vdp[j].Mx * 100, patch.readings[i].Vdp[j].Mx, j)
            # break
        else:
            patch.readings[i].Vdp[j].Mx_err = -99.9
# i = len(patch.readings) - 2
# j = 2
# c, tau, M, err, vs = (
#     jdata.calcColeCole(patch.readings[i].Vdp[j].Vs /
#                        ((patch.readings[i].Vdp[j].Vp)),
#                        patch.window_width,
#                        c_init,
#                        tau_init,
#                        radius))
# End Cole-Cole fitting =============================
stop = clock()
print("timing of cole-cole (sec):", stop)
# print(c, tau, M, err, patch.readings[i].Vdp[j].Mx)
# plt.plot(patch.window_center, -patch.readings[i].Vdp[j].Vs / patch.readings[i].Vdp[j].Vp, 'o-m')
# plt.plot(patch.window_center, vs, 'o-')
# plt.show()
# write to file =====================================
format_jup = '%8s %7s %7s %10s %11.0f %11.0f %14.3f %14.3f %14.3f %14.3f %14.3f %14s %11.0f %14.3f %14.3f %14.3f %13s %12.0f %14.3f %14.3f %14.3f %14.3f %14.3f %14.3f %14.3f %14.3f %14.3f %14.3f %14.3e %14s %11.0f %11s %10.2f %10.2e %7s %8.2f %8.2e %8.2f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n'
out_file = open(outname, "w+")
# write the headers
for i in range(len(patch.headers)):
    if i < (len(patch.headers) - 1):
        out_file.write('%s' % patch.headers[i])
    else:
        out_file.write('%s\n' % patch.headers[i])

# now write the file
for k in range(len(patch.readings)):
    for l in range(len(patch.readings[k].Vdp)):
        out_file.write(format_jup % (str(patch.readings[k].Vdp[l].RDG),
                                     str(patch.readings[k].Vdp[l].dipole),
                                     str(patch.readings[k].Vdp[l].n),
                                     str(patch.readings[k].Vdp[l].Status),
                                     patch.readings[k].Vdp[l].Nom_a,
                                     patch.readings[k].Idp.Tx1,
                                     patch.readings[k].Idp.Tx1East,
                                     patch.readings[k].Idp.Tx1North,
                                     patch.readings[k].Idp.Tx1Elev,
                                     patch.readings[k].Idp.Tx2East,
                                     patch.readings[k].Idp.Tx2North,
                                     str(patch.readings[k].Vdp[l].Rx1File),
                                     patch.readings[k].Vdp[l].Rx1,
                                     patch.readings[k].Vdp[l].Rx1East,
                                     patch.readings[k].Vdp[l].Rx1North,
                                     patch.readings[k].Vdp[l].Rx1Elev,
                                     str(patch.readings[k].Vdp[l].Rx2File),
                                     patch.readings[k].Vdp[l].Rx2,
                                     patch.readings[k].Vdp[l].Rx2East,
                                     patch.readings[k].Vdp[l].Rx2North,
                                     patch.readings[k].Vdp[l].Rx2Elev,
                                     patch.readings[k].Vdp[l].contact,
                                     patch.readings[k].Vdp[l].Vp,
                                     patch.readings[k].Vdp[l].Vp_err,
                                     0.0,
                                     patch.readings[k].Idp.In,
                                     patch.readings[k].Idp.In_err,
                                     patch.readings[k].Vdp[l].Rho,
                                     patch.readings[k].Vdp[l].calcRhoError(patch.readings[k].Idp),
                                     patch.readings[k].Vdp[l].flagRho,
                                     patch.readings[k].Vdp[l].Stack,
                                     patch.readings[k].Vdp[l].TimeBase,
                                     patch.readings[k].Vdp[l].Mx,
                                     patch.readings[k].Vdp[l].Mx_err,
                                     patch.readings[k].Vdp[l].flagMx,
                                     patch.readings[k].Vdp[l].cole_c,
                                     patch.readings[k].Vdp[l].cole_tau,
                                     patch.readings[k].Vdp[l].cole_m,
                                     patch.readings[k].Vdp[l].Vs[0],
                                     patch.readings[k].Vdp[l].Vs[1],
                                     patch.readings[k].Vdp[l].Vs[2],
                                     patch.readings[k].Vdp[l].Vs[3],
                                     patch.readings[k].Vdp[l].Vs[4],
                                     patch.readings[k].Vdp[l].Vs[5],
                                     patch.readings[k].Vdp[l].Vs[6],
                                     patch.readings[k].Vdp[l].Vs[7],
                                     patch.readings[k].Vdp[l].Vs[8],
                                     patch.readings[k].Vdp[l].Vs[9],
                                     patch.readings[k].Vdp[l].Vs[10],
                                     patch.readings[k].Vdp[l].Vs[11],
                                     patch.readings[k].Vdp[l].Vs[12],
                                     patch.readings[k].Vdp[l].Vs[13],
                                     patch.readings[k].Vdp[l].Vs[14],
                                     patch.readings[k].Vdp[l].Vs[15],
                                     patch.readings[k].Vdp[l].Vs[16],
                                     patch.readings[k].Vdp[l].Vs[17],
                                     patch.readings[k].Vdp[l].Vs[18],
                                     patch.readings[k].Vdp[l].Vs[19]))
out_file.close()
print("Done")
