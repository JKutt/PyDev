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


class lokeColumn:
    def __init__(
            self,
            xp,
            yp,
            depths,
            vals):
        self.x = xp
        self.y = yp
        self.elev = depths
        self.val = vals


class ubcColumn:
    def __init__(self,
                 col_cell,
                 xp,
                 yp):
        self.cells = col_cell
        self.x = xp
        self.y = yp
        self.x0 = 0.0
        self.y0 = 0.0
        self.x1 = 0.0
        self.y1 = 0.0
        self.nulled = True


# define the file required for import
fileName = "/Users/juan/Documents/testData/Fremont_Culled_Loke_rotated.xyz"
outName = "/Users/juan/Documents/testData/Fremont_Culled_Loke-trim.con"
outMesh = "/Users/juan/Documents/testData/Fremont_Culled_Loke-trim.msh"
lines = 0
splarray = " "

# set variables for
trans_point_X = 753756.00
trans_point_Y = 4163372.0
depth_cut = 200.0
theta = 22.0 * np.pi / 180
maxE = 701945.0
maxN = 7395445.0
minE = 600345.0
minN = 6389445.0
null_flag = 1e-8
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
x = np.zeros(lines - 1)                 # number of x cells
y = np.zeros(lines - 1)                 # number of y cells
z = np.zeros(lines - 1)                 # number of z cells
rho = np.zeros(lines - 1)               # array of values
rot = np.array([[np.cos(theta), -np.sin(theta)],
               [np.sin(theta), np.cos(theta)]])
print(rot)
currRdg = 0
cnt = 0
dpt = np.zeros(((lines - 1), 3))
for i, line in enumerate(text_file):
        Varinfo = line.split()
        if i > 0:
            x[cnt] = float(Varinfo[0])
            y[cnt] = float(Varinfo[1])
            z[cnt] = float(Varinfo[2])
            rho[cnt] = float(Varinfo[3])
            dpt[cnt, 0] = x[cnt]
            dpt[cnt, 1] = y[cnt]
            dpt[cnt, 2] = z[cnt]
            cnt += 1
x_sort = np.sort(x)
y_sort = np.sort(y)
# y.sort()
# z.sort()
d_x = []                            # x cell spacing
d_y = []                            # y cell spacing
d_z = []
z_cell = []
dist_x = 0
dist_y = 0
for i in range(x.size - 1):
    distx = x_sort[i + 1] - x_sort[i]
    if distx > 0:
        d_x.append(distx)
        dist_x = distx
    disty = y_sort[i + 1] - y_sort[i]
    if disty > 0:
        d_y.append(disty)
        dist_y = disty
    distz = z[i + 1] - z[i]
    if distz > 0:
        d_z.append(distz)

for i in range(x.size):
    if ((x[i] == x[0]) & (y[i] == y[0])):
        z_cell.append(z[i])

print np.min(x), np.max(x)
print np.min(y), np.max(y)
print np.min(z), np.max(z)
# print dpt
num_x_cells = len(d_x) + 1
num_y_cells = len(d_y) + 1

# I feel like I have to build the x and y stations
curr_stnx = 0.0
cntr = 0
x_stns = np.zeros(num_x_cells)
y_stns = np.zeros(num_y_cells)
x_stns[0] = np.min(x)
y_stns[0] = np.min(y)
for i in range(1, num_x_cells):
    x_stns[i] = x_stns[0] + sum(d_x[0:i])
for i in range(1, num_y_cells):
    y_stns[i] = y_stns[0] + sum(d_y[0:i])

#####################################################
# now loop through and assign all columns to the given x and y
col_cell = []
for r in range(y_stns.size):
    for n in range(x_stns.size):
        idx_x = np.where(x == x_stns[n])[0]
        idx_y = np.where(y == y_stns[r])[0]
        c_z = []
        c_v = []
        for i in range(idx_x.size):
            if y[idx_x[i]] == y_stns[r]:
                c_z.append(z[idx_x[i]])
                c_v.append(rho[idx_x[i]])
        col = lokeColumn(x_stns[n], y_stns[r], c_z, c_v)
        col_cell.append(col)

#####################################################
# start creating the z cells given all new info
elev_ = np.zeros(len(col_cell))
for i in range(len(col_cell)):
    elev_[i] = col_cell[i].elev[0]

min_elev = np.min(elev_)
max_elev = np.max(elev_)
zsize = np.zeros(len(col_cell[0].elev) - 1)
for i in range(len(col_cell[0].elev) - 1):
    zsize[i] = col_cell[0].elev[i] - col_cell[0].elev[i + 1]

first_cell_min = min_elev + (zsize[0] / 2.0)
first_cell_max = max_elev + (zsize[0] / 2.0)
num_extra_cells = (first_cell_max - first_cell_min) / zsize[0]
total_num_zcell = num_extra_cells + len(col_cell[0].elev)
print first_cell_max, first_cell_min, total_num_zcell

cell_div = []
cntt = 0
for k in range(int(total_num_zcell)):
    if k < num_extra_cells:
        cell_div.append(zsize[0])
    else:
        cell_div.append(zsize[cntt])
        cntt += 1
cell_div.append(zsize[cntt - 1])
# print cell_div
levels = np.zeros(len(cell_div))
levels[0] = first_cell_max
for i in range(1, len(cell_div)):
    levels[i] = levels[0] - sum(cell_div[0:i])

##################################################
# create the ubc mesh cells
ubcColCells = []
for n in range(len(col_cell)):
    depth_filter = 0.0                         # depth filter for each location
    first_elev = False                         # flag model topo
    ubc_ = []
    xp = (
        col_cell[n].x * np.cos(theta)) + (
        col_cell[n].y * np.sin(theta))
    yp = (
        -col_cell[n].x * np.sin(theta)) + (
        col_cell[n].y * np.cos(theta))
    for j in range(levels.size - 1):
        check = False                          # flag for NULL cell
        for r in range(len(col_cell[n].elev)):
            if (levels[j] > col_cell[n].elev[r] > levels[j + 1]) and (col_cell[n].elev[r] > depth_filter):
                # print "found cell: ", col_cell[n].elev[r]
                # print " @: ", levels[j], levels[j + 1]
                ubc_.append(col_cell[n].val[r])
                check = True
                break
        if check:
            if first_elev:
                pass
            else:
                first_elev = True
                depth_filter = levels[j] - depth_cut
                # print "max depth for column: ", depth_filter
        else:
            # print " NULL @: ", levels[j], levels[j + 1]
            ubc_.append(null_flag)
    ubcColCells.append(ubcColumn(ubc_, xp, yp))        # assigns the ubc column

#################################################
# rotate the data

xp = np.zeros(len(ubcColCells))
yp = np.zeros(len(ubcColCells))
maxX = ubcColCells[0].x
minX = ubcColCells[0].x
maxY = ubcColCells[0].y
minY = ubcColCells[0].y
# maxZ = ubcColCells[0].z
# minZ = ubcColCells[0].z
for j in range(len(ubcColCells)):
    xp[j] = ubcColCells[j].x
    yp[j] = ubcColCells[j].y
    if maxX < ubcColCells[j].x:
        maxX = ubcColCells[j].x
    if minX > ubcColCells[j].x:
        minX = ubcColCells[j].x
    if maxY < ubcColCells[j].y:
        maxY = ubcColCells[j].y
    if minY > ubcColCells[j].y:
        minY = ubcColCells[j].y
    # if maxZ < ubcColCells[j].z:
    #     maxZ = ubcColCells[j].z
    # if minz > ubcColCells[j].z:
        # minz = ubcColCells[j].z
#################################################
# create new mesh and model to new rotated data
half_xcell = dist_x / 4.0
half_ycell = dist_y / 4.0
numZcell = len(ubcColCells[0].cells)
numXcell_new = np.ceil(
    ((maxX + half_xcell) - (minX - half_xcell)) / (dist_x / 2.0))
numYcell_new = np.ceil(
    ((maxY + half_ycell) - (minY - half_ycell)) / (dist_y / 2.0))
header1 = str(numXcell_new) + ' ' + str(numYcell_new) + ' ' + str(numZcell)

new_cells = []
x_length = numXcell_new * (dist_x / 2.0)
x_new = np.arange(0, x_length, (dist_x / 2.0))
y_length = numYcell_new * (dist_y / 2.0)
y_new = np.arange(0, y_length, (dist_y / 2.0))
header3 = str(numXcell_new) + '*' + str(dist_x / 2.0)
header4 = str(numYcell_new) + '*' + str(dist_y / 2.0)
# lets create a null column
nulcol = ubcColCells[0].cells
for i in range(len(nulcol)):
    nulcol[i] = null_flag
# create the new models cells
for i in range(y_new.size - 1):
    for j in range(x_new.size - 1):
        mcell = ubcColumn(nulcol, 0.0, 0.0)
        mcell.x0 = x_new[j]
        mcell.x1 = x_new[j + 1]
        mcell.y0 = y_new[i]
        mcell.y1 = y_new[i + 1]
        new_cells.append(mcell)

# search through and assign data
for idx in range(len(new_cells)):
    for idx1 in range(len(ubcColCells)):
        if ((new_cells[idx].x0 < ubcColCells[idx1].x < new_cells[idx].x1) and (new_cells[idx].y0 < ubcColCells[idx1].y < new_cells[idx].y1)):
            if new_cells[idx].nulled:
                new_cells[idx].cells = ubcColCells[idx1].cells
                new_cells[idx].nulled = False
                # print new_cells[idx].cells
            else:
                for k in range(len(new_cells[idx].cells)):
                    new_cells[idx].cells[k] = (
                        new_cells[idx].cells[k] +
                        ubcColCells[idx1].cells[k]) / 2.0


# likely a waiste but just to get it done; translate back to proper coords.
min_e = new_cells[0].x0 + trans_point_X
min_n = new_cells[0].y0 + trans_point_Y
max_z = levels[0]
for idx in range(len(new_cells)):
    new_cells[idx].x0 = new_cells[idx].x0 + trans_point_X
    new_cells[idx].x1 = new_cells[idx].x1 + trans_point_X
    new_cells[idx].y0 = new_cells[idx].y0 + trans_point_Y
    new_cells[idx].y1 = new_cells[idx].y1 + trans_point_Y
    if min_e > new_cells[idx].x0:
        min_e = new_cells[idx].x0
    if min_n > new_cells[idx].y0:
        min_n = new_cells[idx].y0

header2 = str(min_e) + ' ' + str(min_n) + ' ' + str(max_z)
#################################################
# and done! now write it to file
out_file = open(outName, "w")
# now the formatting has been completed
for idx in range(len(new_cells)):
    for k in range(len(new_cells[idx].cells)):
        try:
            out_file.write("%0.5e\n" % 1.0 / new_cells[idx].cells[k])
        except:
            out_file.write("%0.5e\n" % 1e-8)

out_file.close()

header5 = ''
for i in range(len(cell_div)):
    if i == 0:
        header5 = str(cell_div[i])
    else:
        header5 = header5 + " " + str(cell_div[i])

# write meshg file
out_msh_file = open(outMesh, "w")
out_msh_file.write("%s\n" % header1)
out_msh_file.write("%s\n" % header2)
out_msh_file.write("%s\n" % header3)
out_msh_file.write("%s\n" % header4)
out_msh_file.write("%s\n" % header5)
out_msh_file.close()
# print len(col)
plt.plot(xp, yp, 'o')
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(x, y, z, c='r')
plt.show()
