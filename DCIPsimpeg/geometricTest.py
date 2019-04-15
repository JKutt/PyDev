import matplotlib.pyplot as plt
import numpy as np
import DCIPtools as DCIP
import matplotlib.colors as colors


def calcHomoHalfSpaceVoltage(background, r1, r2):
    In = 1                                           # normalized current
    v1 = (In * background) / (2 * np.pi * r1)        # voltage at point 1
    v2 = -1 * (In * background) / (2 * np.pi * r2)      # voltage at point 2
    v = v1 + v2                                      # theoretical voltage
    return v


# ============================ User area ==============================
# set variables
fname = "/Users/juan/Documents/testData/B2_15_125-Loc.DAT"
outname = "/Users/juan/Documents/testData/B2_15_125-Loc-mod.DAT"

# tx1 = [388765.000, 3279595.000, 28.000]
# tx2 = [388755.000, 3279553.000, 31.000]
# rx1 = [388710.000, 3279591.000, 26.000]
# rx2 = [388713.000, 3279607.000, 26.000]

# tx1 = [403499.407, 4418560.000, 1927.]
# tx2 = [411852.000, 4419621.000, 2000.]
# rx1 = [403861.830, 4421202.177, 5084.]
# rx2 = [403875.037, 4421300.810, 2147.]

# tx1 = [389414., 3279056., 20]
# tx2 = [389428., 3279114., 20]
# rx1 = [389378., 3279097., 20]
# rx2 = [389345., 3279002., 20]

# tx1 = [388755.000, 3279553.000, 31.000]
# tx2 = [388765.000, 3279595.000, 28.000]
# rx1 = [388673.000, 3279494.000, 28.000]
# rx2 = [388571.000, 3279507.000, 28.000]

# tx1 = [388765.000, 3279595.000, 28.000]
# tx2 = [388755.000, 3279553.000, 31.000]
# rx1 = [388593.000, 3279577.000, 28.000]
# rx2 = [388571.000, 3279507.000, 28.000]

# tx1 = [389428.000, 3279114.000, 20.000]
# tx2 = [389414.000, 3279056.000, 20.000]
# rx1 = [389390.000, 3279120.000, 18.000]
# rx2 = [389396.000, 3279142.000, 18.000]

tx1 = [402498.136, 4422553.847, 1954.000]
tx2 = [402551.309, 4418553.113, 1843.000]
# rx1 = [407872.933, 4416797.805, 2248.700]
# rx2 = [407624.695, 4416704.227, 2232.320]

rx1 = [402626.000, 4419400.000, 1850.000]
rx2 = [402128.000, 4419102.000, 1821.000]

# points = np.vstack((tx2, tx1, rx1, rx2))
points = np.vstack((tx1, tx2, rx1, rx2))
minx = np.min(points[:, 0]) - 100
miny = np.min(points[:, 1]) - 100
maxx = np.max(points[:, 0]) + 100
maxy = np.max(points[:, 1]) + 100
print("min X: {0}".format(minx))
print("max X: {0}".format(maxx))
print("min Y: {0}".format(miny))
print("max Y: {0}".format(maxy))
x = np.linspace(minx, maxx, 25)
y = np.linspace(miny, maxy, 25)
X, Y = np.meshgrid(x, y)
tx1x = np.zeros_like(X)
tx1y = np.zeros_like(Y)
tx1x = tx1x + tx1[0]
tx1y = tx1y + tx1[1]
tx2x = np.zeros_like(X)
tx2y = np.zeros_like(Y)
tx2x = tx2x + tx2[0]
tx2y = tx2y + tx2[1]

r1_ = ((X - tx1x)**2 + (Y - tx1y)**2)**0.5
r2_ = ((X - tx2x)**2 + (Y - tx2y)**2)**0.5
In = 1
background = 5000                                 # normalized current
v1 = (In * background) / (2 * np.pi * r1_)        # voltage at point 1
v2 = -1 * (In * background) / (2 * np.pi * r2_)      # voltage at point 2
v = (v1 + v2) * 1000
plt.contourf(X, Y, v, cmap='magma')
plt.colorbar()
# for idx in range(x.size):

r1 = ((rx1[0] - tx1[0])**2 +
      (rx1[1] - tx1[1])**2 +
      (rx1[2] - tx1[2])**2)**0.5
r2 = ((rx2[0] - tx1[0])**2 +
      (rx2[1] - tx1[1])**2 +
      (rx2[2] - tx1[2])**2)**0.5
r3 = ((rx1[0] - tx2[0])**2 +
      (rx1[1] - tx2[1])**2 +
      (rx1[2] - tx2[2])**2)**0.5
r4 = ((rx2[0] - tx2[0])**2 +
      (rx2[1] - tx2[1])**2 +
      (rx2[2] - tx2[2])**2)**0.5

gf = 1 / ((1 / r1 - 1 / r2) - (1 / r3 - 1 / r4))

v1_ = calcHomoHalfSpaceVoltage(background, r1, r3)
v2_ = calcHomoHalfSpaceVoltage(background, r2, r4)
V_ = v1_ - v2_
print("Voltage at node 1: {0} & node 2: {1}".format(v1_, v2_))
print("Dipole theretical Voltage: {0}".format(V_))
print("G: {0}".format(2 * np.pi * gf))
print("r1: {0} r2: {1} r3: {2} r4: {3}".format(r1, r2, r3, r4))
# if r1 < r3:
#     if gf < 0:
#         if r1 > r2:
#             temp_rx1 = rx1
#             temp_rx2 = rx2
#             rx1 = temp_rx2
#             rx2 = temp_rx1
# else:
#     if gf < 0:
#         temp_tx1 = tx1
#         temp_tx2 = tx2
#         tx1 = temp_tx2
#         tx2 = temp_tx1

# r1 = ((rx1[0] - tx1[0])**2 +
#       (rx1[1] - tx1[1])**2 +
#       (rx1[2] - tx1[2])**2)**0.5
# r2 = ((rx2[0] - tx1[0])**2 +
#       (rx2[1] - tx1[1])**2 +
#       (rx2[2] - tx1[2])**2)**0.5
# r3 = ((rx1[0] - tx2[0])**2 +
#       (rx1[1] - tx2[1])**2 +
#       (rx1[2] - tx2[2])**2)**0.5
# r4 = ((rx2[0] - tx2[0])**2 +
#       (rx2[1] - tx2[1])**2 +
#       (rx2[2] - tx2[2])**2)**0.5

# gf = 1 / ((1 / r1 - 1 / r2) - (1 / r3 - 1 / r4))

# print("G: {0}".format(2 * np.pi * gf))
# print("r1: {0} r2: {1} r3: {2} r4: {3}".format(r1, r2, r3, r4))

plt.plot(tx1[0], tx1[1], 'r*', label='tx1')
plt.plot(tx2[0], tx2[1], 'k*', label='tx2')
plt.plot(rx1[0], rx1[1], 'bo', label='rx1')
plt.plot(rx2[0], rx2[1], 'go', label='rx2')
plt.title("Gf: {:f}".format(2 * np.pi * gf))
plt.legend()
plt.axes().set_aspect('equal', 'datalim')
plt.show()
