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
fileName = "C:/Users/johnk/Projects/Karst/loke/B10-DC-iter13.xyz"
outName = "C:/Users/johnk/Projects/Karst/loke/B10-DC-iter13-trim.con"
outMesh = "C:/Users/johnk/Projects/Karst/loke/B10-DC-iter13-trim.msh"
lines = 0
splarray = " "

# set variables for
trans_point_X = 753756.00
trans_point_Y = 4163372.0
depth_cut = 200.0
theta = 0
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

    