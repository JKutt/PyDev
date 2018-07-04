import matplotlib.pyplot as plt
import numpy as np
import JDataObject as Jdata

# define the file required for import
fileName = "/Users/juan/Documents/testData/L6050-err1.DAT"
fileName2 = "/Users/juan/Documents/testData/L6050-err2.DAT"
lines = 0
splarray = " "

# designate format
outFormat = [7, 8, 10, 11, 14, 14, 14, 14, 14, 14, 14, 14, 14, 13, 14,
    14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 11, 11, 10, 10, 7, 9,
     9, 9, 9, 9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9]

# open the file 1 ------------------------------------------------
text_file = open(fileName, "r")

# determin how many lines in the file
while text_file.readline():
        lines += 1
print(lines)
text_file.close()

# initiate a patch
patch = Jdata.Jpatch(6051)
patches = []
patches.append(patch)
# read header information
text_file = open(fileName, "r")
# initiate reading control variables
currRdg = 0
for i, line in enumerate(text_file):
    if i == 4:
        Varinfo = line.split()
        header4 = line
        # print(Varinfo)
    elif i == 0:
            header1 = line
    elif i == 1:
            header2 = line
    elif i == 2:
            header3 = line
    elif i > 3:
        # try:
                datatxt = line.split()
                # do some Jdatamanagment stuff
                varFields = Jdata.Jreadtxtline(Varinfo, datatxt)
                # verify if line is a new reading
                if varFields.RDG == currRdg:
                    # add the dipoles
                    # Idp = Jdata.JinDipole(varFields)
                    Vpdp = Jdata.JvoltDipole(varFields)
                    Rdg.addVoltageDipole(Vpdp)
                else:
                    # create a reading
                    Rdg = Jdata.Jreading(varFields.RDG)
                    Idp = Jdata.JinDipole(varFields)
                    Vpdp = Jdata.JvoltDipole(varFields)
                    Rdg.addVoltageDipole(Vpdp)
                    Rdg.addInDipole(Idp)
                    # add reading to the patch
                    patches[0].addreading(Rdg)
                    currRdg = varFields.RDG
        # except:
                # pass

text_file.close()

# Open second file --------------------------------------
text_file = open(fileName, "r")
# determin how many lines in the file
while text_file.readline():
        lines += 1
print(lines)
text_file.close()
# initiate a second patch
patch2 = Jdata.Jpatch(6052)
patches.append(patch2)
# read header information
text_file = open(fileName2, "r")
# initiate reading control variables
currRdg = 0
for i, line in enumerate(text_file):
    if i == 4:
        Varinfo = line.split()
        header4 = line
        # print(Varinfo)
    elif i == 0:
            header1 = line
    elif i == 1:
            header2 = line
    elif i == 2:
            header3 = line
    elif i > 3:
        # try:
                datatxt = line.split()
                # do some Jdatamanagment stuff
                varFields = Jdata.Jreadtxtline(Varinfo, datatxt)
                # verify if line is a new reading
                if varFields.RDG == currRdg:
                    # add the dipoles
                    # Idp = Jdata.JinDipole(varFields)
                    Vpdp = Jdata.JvoltDipole(varFields)
                    Rdg.addVoltageDipole(Vpdp)
                else:
                    # create a reading
                    Rdg = Jdata.Jreading(varFields.RDG)
                    Idp = Jdata.JinDipole(varFields)
                    Vpdp = Jdata.JvoltDipole(varFields)
                    Rdg.addVoltageDipole(Vpdp)
                    Rdg.addInDipole(Idp)
                    # add reading to the patch
                    patches[1].addreading(Rdg)
                    currRdg = varFields.RDG
        # except:
                # pass
text_file.close()

# Plotting ---------------------------------------------
# print(Varinfo[1])
print("chosendipole info:")
print(len(patches[1].readings[1].Vdp))
num_wins = patches[1].readings[1].Vdp[0].Vs.size
# Window times
timeFrom = [2040, 2060, 2080, 2120, 2160, 2200,
        2240, 2320, 2400,
        2480, 2560, 2640,
        2720, 2800, 2960,
        3120, 3280, 3440,
        3600, 3760]
timeTo = [2060, 2080, 2120, 2160, 2200, 2240,
        2320,2400,2480,2560,2640,2720
        ,2800,2960,3120,3280,3440
        ,3600,3760,3920]
t_ = np.zeros(num_wins)
for i in range(num_wins):
    t_[i] = (timeFrom[i] + timeTo[i]) / 2.0

f, plts = plt.subplots(2, 3)
plts[0, 0].plot(t_, -patches[0].readings[0].Vdp[1].Vs, 'b-o')
plts[0, 0].plot(t_, -patches[1].readings[0].Vdp[1].Vs, 'm-o')
plts[0, 0].set_title("Mem164 blue = Weighted & Magenta = Averaged")
plts[0, 1].plot(t_, -patches[0].readings[0].Vdp[2].Vs, 'b-o')
plts[0, 1].plot(t_, -patches[1].readings[0].Vdp[2].Vs, 'm-o')
plts[0, 2].plot(t_, -patches[0].readings[0].Vdp[3].Vs, 'b-o')
plts[0, 2].plot(t_, -patches[1].readings[0].Vdp[3].Vs, 'm-o')
plts[1, 0].plot(t_, patches[0].readings[0].Vdp[4].Vs, 'b-o')
plts[1, 0].plot(t_, patches[1].readings[0].Vdp[4].Vs, 'm-o')
plts[1, 1].plot(t_, -patches[0].readings[0].Vdp[5].Vs, 'b-o')
plts[1, 1].plot(t_, -patches[1].readings[0].Vdp[5].Vs, 'm-o')
plts[1, 2].plot(t_, -patches[0].readings[0].Vdp[6].Vs, 'b-o')
plts[1, 2].plot(t_, -patches[1].readings[0].Vdp[6].Vs, 'm-o')
plt.show()
# print(patches[0].readings[10].Vdp[2].Rx1)
# print(patches[0].readings[10].Idp.Tx1)
# print(patches[0].readings[10].Idp.In)