import matplotlib.pyplot as plt
import numpy as np
import JDataObject as Jdata

# define the file required for import
fileName = "C:/Users/johnk/Projects/nothernStar/18060NSTA_8780/18060NSTA_8800.ip"

sm24_start_times = []
sm24_stop_times = []
datalines = []
lines = 0
text_file = open(fileName, "r")

# determin how many lines in the file
while text_file.readline():
        lines += 1
text_file.close()

# initiate a patch
patch = Jdata.Jpatch()
# read header information
text_file = open(fileName, "r")
# initiate reading control variables
currRdg = 0
for i, line in enumerate(text_file):
    if i == 2:
        Varinfo = line.split()
        header4 = line
        # print(Varinfo)
    elif i == 3:
        extract = line.split('=')
        extract2 = extract[1].split(',')
        sm24_start_times = extract2
    elif i == 4:
        extract = line.split('=')
        extract2 = extract[1].split(',')
        sm24_stop_times = extract2
    elif i > 11:
        try:
            datatxt = line.split()
            # # do some Jdatamanagment stuff
            # # print(datatxt)
            varFields = Jdata.Jreadtxtline(Varinfo, datatxt)
            datalines.append(varFields)
            # # verify if line is a new reading
            # if varFields.RDG == currRdg:
            #     # add the dipoles
            #     # Idp = Jdata.JinDipole(varFields)
            #     Vpdp = JvoltDipole(varFields)
            #     Rdg.addVoltageDipole(Vpdp)
            # else:
            #     # create a reading
            #     Rdg = Jreading(varFields.RDG)
            #     Idp = JinDipole(varFields)
            #     Vpdp = JvoltDipole(varFields)
            #     Rdg.addVoltageDipole(Vpdp)
            #     Rdg.addInDipole(Idp)
            #     # add reading to the patch
            #     patch.addreading(Rdg)
            #     currRdg = varFields.RDG
        except:
            pass
# print(len(datalines))
# print(datalines[0].RDNG)
sm24_start_times = np.asarray(sm24_start_times).astype(float)
sm24_stop_times = np.asarray(sm24_stop_times).astype(float)
# print(sm24_start_times)
# print(sm24_stop_times)
mid_points = (sm24_start_times + sm24_stop_times) / 2.
# rdg = 11
# for index in range(len(datalines)):
#     if rdg == int(datalines[index].RDNG):
#         decay = [datalines[index].CH1,
#                  datalines[index].CH2,
#                  datalines[index].CH3,
#                  datalines[index].CH4,
#                  datalines[index].CH5,
#                  datalines[index].CH6,
#                  datalines[index].CH7,
#                  datalines[index].CH8,
#                  datalines[index].CH9,
#                  datalines[index].CH10,
#                  datalines[index].CH11]
#         plt.plot(mid_points, decay, 'r-o')
# plt.show()
index = 17
decay = [datalines[index].CH1,
         datalines[index].CH2,
         datalines[index].CH3,
         datalines[index].CH4,
         datalines[index].CH5,
         datalines[index].CH6,
         datalines[index].CH7,
         datalines[index].CH8,
         datalines[index].CH9,
         datalines[index].CH10,
         datalines[index].CH11]
plt.plot(mid_points, decay, 'r-o')
decay = [datalines[index + 4].CH1,
         datalines[index + 4].CH2,
         datalines[index + 4].CH3,
         datalines[index + 4].CH4,
         datalines[index + 4].CH5,
         datalines[index + 4].CH6,
         datalines[index + 4].CH7,
         datalines[index + 4].CH8,
         datalines[index + 4].CH9,
         datalines[index + 4].CH10,
         datalines[index + 4].CH11]
plt.plot(mid_points, decay, 'r-o')
decay = [datalines[index + 8].CH1,
         datalines[index + 8].CH2,
         datalines[index + 8].CH3,
         datalines[index + 8].CH4,
         datalines[index + 8].CH5,
         datalines[index + 8].CH6,
         datalines[index + 8].CH7,
         datalines[index + 8].CH8,
         datalines[index + 8].CH9,
         datalines[index + 8].CH10,
         datalines[index + 8].CH11]
plt.plot(mid_points, decay, 'r-o')
print(datalines[index].RDNG)
print(datalines[index].C1E)
print(datalines[index].P1E)
print(datalines[index].P2E)
print(datalines[index + 4].RDNG)
print(datalines[index + 4].C1E)
print(datalines[index + 4].P1E)
print(datalines[index + 4].P2E)
print(datalines[index + 8].RDNG)
print(datalines[index + 8].C1E)
print(datalines[index + 8].P1E)
print(datalines[index + 8].P2E)

# define the file required for import
fileName = "C:/Users/johnk/Projects/nothernStar/18060NSTA_8780/L8900-test.DAT"

patch = Jdata.loadDias(fileName)   # Create the patch from data file
rdg = 19                           # Source to plot
dp = 28
# i = 4
dpnum = ['SM24','SM24','SM24', 'DIAS32 - 2', 'DIAS32 - 1','DIAS32 - 3', 'DIAS32 - 4', 'DIAS32 - Ens']
print(len(patch.readings))
print(patch.readings[rdg].Idp.Tx1)
print(patch.readings[rdg].Vdp[dp].Rx1)
print(patch.readings[rdg].Vdp[dp].Rx2)
plt.plot(patch.window_center,
         patch.readings[rdg].Vdp[dp].Vs / (patch.readings[rdg].Vdp[dp].Vp / 1000.), 'k-o')

# define the file required for import
fileName2 = "C:/Users/johnk/Projects/nothernStar/18060NSTA_8780/L8900-test-1.DAT"

patch2 = Jdata.loadDias(fileName2)   # Create the patch from data file
rdg = 19                           # Source to plot
dp = 28
# i = 4
print(len(patch2.readings))
print(patch2.readings[rdg].Idp.Tx1)
print(patch2.readings[rdg].Vdp[dp].Rx1)
print(patch2.readings[rdg].Vdp[dp].Rx2)
plt.plot(patch2.window_center,
         patch2.readings[rdg].Vdp[dp].Vs / (patch2.readings[rdg].Vdp[dp].Vp / 1000.), 'm-o')

# define the file required for import
fileName3 = "C:/Users/johnk/Projects/nothernStar/18060NSTA_8780/L8900-test-3.DAT"

patch3 = Jdata.loadDias(fileName3)   # Create the patch from data file
rdg = 19                             # Source to plot
dp = 28
# i = 4
print(len(patch3.readings))
print(patch3.readings[rdg].Idp.Tx1)
print(patch3.readings[rdg].Vdp[dp].Rx1)
print(patch3.readings[rdg].Vdp[dp].Rx2)
plt.plot(patch3.window_center,
         patch3.readings[rdg].Vdp[dp].Vs / (patch3.readings[rdg].Vdp[dp].Vp / 1000.), 'c-o')

# define the file required for import
fileName4 = "C:/Users/johnk/Projects/nothernStar/18060NSTA_8780/L8900-test-4.DAT"

patch4 = Jdata.loadDias(fileName4)   # Create the patch from data file
rdg = 19                           # Source to plot
dp = 28
# i = 4
print(len(patch4.readings))
print(patch4.readings[rdg].Idp.Tx1)
print(patch4.readings[rdg].Vdp[dp].Rx1)
print(patch4.readings[rdg].Vdp[dp].Rx2)
plt.plot(patch4.window_center,
         patch4.readings[rdg].Vdp[dp].Vs / (patch4.readings[rdg].Vdp[dp].Vp / 1000.), 'g-o')

# define the file required for import
fileName5 = "C:/Users/johnk/Projects/nothernStar/18060NSTA_8780/L8900-test-4ens.DAT"

patch5 = Jdata.loadDias(fileName5)   # Create the patch from data file
rdg = 19                           # Source to plot
dp = 28
# i = 4
print(len(patch5.readings))
print(patch5.readings[rdg].Idp.Tx1)
print(patch5.readings[rdg].Vdp[dp].Rx1)
print(patch5.readings[rdg].Vdp[dp].Rx2)
plt.plot(patch5.window_center,
         patch5.readings[rdg].Vdp[dp].Vs / (patch5.readings[rdg].Vdp[dp].Vp / 1000.), 'b-o')

plt.legend(dpnum)
plt.title("Tx 3600 - 3800 Rx 5700 - 5800")
plt.xlabel("time (ms)")
plt.ylabel("Voltage (mV/V)")
plt.show()
