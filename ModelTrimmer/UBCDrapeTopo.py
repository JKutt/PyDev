import matplotlib.pyplot as plt
import numpy as np


class modelCell:
    # class for UBC inversion model cell
    def __inti__(self):
        self.x0 = 0
        self.x1 = 0
        self.y0 = 0
        self.y1 = 0
        self.z0 = 0
        self.z1 = 0
        self.val = 0


# define the file required for import
# fileName0 = "/home/juan/Documents/DevProjects/testData/Resampling/DEM.txt" 
fileName = "/home/juan/Documents/DevProjects/testData/Resampling/MergedMesh_Fine.txt"
fileName2 = "/home/juan/Documents/DevProjects/testData/Resampling/ResampledCoarseInversion_Iteration12.con"
outName = "/home/juan/Documents/DevProjects/testData/Resampling/dcinv3d-trim.con"
lines = 0
splarray = " "

# set variables for
# maxE = 701945.0
# maxN = 7395445.0
# minE = 600345.0
# minN = 6389445.0
# max_depth = 100.0
null_flag = 1e-8
# open file 0 ====================================================
# csv = open(fileName0)
# easting = []
# northing = []
# dem = []
# for dline in csv:
#     splitt = dline.split(" ")
#     if len(splitt) > 1:
#         easting.append(splitt[0])
#         northing.append(splitt[1])
#         dem.append(splitt[2])
# csv.close()

# open the file 1 ------------------------------------------------
text_file = open(fileName, "r")

# determin how many lines in the file
while text_file.readline():
        lines += 1
print(lines)
text_file.close()

# read header information
text_file = open(fileName, "r")
# initiate reading control variables
xcell = 0                              # number of x cells
ycell = 0                              # number of y cells
zcell = 0                              # number of z cells
x0 = 0.0                               # initiate corner
y0 = 0.0                               # initiate corner
z0 = 0.0                               # initiate corner
x_format = []                          # inititate formats
y_format = []                          # inititate formats
z_format = []                          # inititate formats
currRdg = 0
for i, line in enumerate(text_file):
        Varinfo = line.split()
        if i == 0:
            xcell = float(Varinfo[0])
            ycell = float(Varinfo[1])
            zcell = float(Varinfo[2])
        elif i == 1:
            x0 = float(Varinfo[0])
            y0 = float(Varinfo[1])
            z0 = float(Varinfo[2])
        elif i == 2:
            x_format = Varinfo
        elif i == 3:
            y_format = Varinfo
        elif i == 4:
            z_format = Varinfo

# open file 2 ------------------------------------------------
text_file = open(fileName2, "r")
lines = 0
# determin how many lines in the file
while text_file.readline():
        lines += 1
print(lines)
text_file.close()

# read header information
text_file = open(fileName2, "r")
# initiate reading control variables
value = np.zeros(lines)                              # number of x cells
currRdg = 0
for i, line in enumerate(text_file):
    value[i] = float(line)
text_file.close()

# now that we got all the info lets try constructing the model
# start with finding x - y - z cell nodes
x_tell = []
for i in range(len(x_format)):
        split = x_format[i].split('*')
        if len(split) == 1:
            x_tell.append(float(x_format[i]))
        elif len(split) > 1:
            number_of_cell = int(split[0])
            for j in range(number_of_cell):
                x_tell.append(float(split[1]))
y_tell = []
for i in range(len(y_format)):
        split = y_format[i].split('*')
        if len(split) == 1:
            y_tell.append(float(y_format[i]))
        elif len(split) > 1:
            number_of_cell = int(split[0])
            for j in range(number_of_cell):
                y_tell.append(float(split[1]))
z_tell = []
for i in range(len(z_format)):
        split = z_format[i].split('*')
        if len(split) == 1:
            z_tell.append(float(z_format[i]))
        elif len(split) > 1:
            number_of_cell = int(split[0])
            for j in range(number_of_cell):
                z_tell.append(float(split[1]))
# now lets build the first row of cells #########
cell = []
x_ = np.zeros(len(x_tell))
y_ = np.zeros(len(y_tell))
z_ = np.zeros(len(z_tell))
cnt = 0
# I think can be removed
for i in range(len(x_tell)):
    x_[i] = x_tell[i]
for i in range(len(y_tell)):
    y_[i] = y_tell[i]
for i in range(len(z_tell)):
    z_[i] = z_tell[i]

########################################
cnt = 0
for k in range(len(x_)):
    for j in range(len(y_)):
        first_lev = False
        depth_ref = 0
        hold_first_value = 0
        cntr = cnt
        for i in range(z_.size):
            if (value[cntr] != 1e-8) and (first_lev == False):
                depth_ref = i
                first_lev = True
                hold_first_value = value[cntr]
                break 
            cntr += 1           
        print("first level is: ", depth_ref)
        for i in range(z_.size):
            mcell = modelCell()
            mcell.x0 = x0 + np.sum(x_[0:j])
            mcell.y0 = y0 + np.sum(y_[0:k])
            mcell.z0 = z0 - np.sum(z_[0:i])
            mcell.x1 = mcell.x0 + x_[j]
            mcell.y1 = mcell.y0 + y_[k]
            mcell.z1 = mcell.z0 - z_[i]
            if i < depth_ref:
                mcell.val = hold_first_value
            else:
                mcell.val = value[cnt]
            # print mcell.val
            cell.append(mcell)
            cnt += 1

# need to print to file now
out_file = open(outName, "w")
# now the formatting has been completed
for idx in range(len(cell)):
    out_file.write("%0.5e\n" % cell[idx].val)

out_file.close()
# SUCCESS !!!!!!!!!
