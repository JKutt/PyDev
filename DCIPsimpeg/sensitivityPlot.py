import numpy as np
import DCIPtools as DCIP
import matplotlib.pyplot as plt


# ============================ User area ==============================
# ver 0.1
# set variables
# exapanding dipoles -------------------------------
dx = 12.5
dx_rx = 50.0
xb = -6000.0
yb = 0.0
zb = 1000.0
xa = -1400.0
ya = 0.0
za = 1000.0
xm_0 = np.arange(-1400.0, -750.0, 50)
xm_1 = np.arange(-700, 0, 100)
xm_2 = np.arange(0, 800, 200)
xm_3 = np.arange(600, 1800, 400)
# xm_4 = np.arange(400, 800, 200)
# xm_5 = np.arange(1000, 1800, 400)
xm = np.hstack((xm_0, xm_1, xm_2, xm_3))
# # xm = np.arange(-1000.0, 1050.0, dx_rx)
ym = 0.0
zm = 1000.0
xn = np.zeros(xm.size)
xn[:-1] = xm[1:]
xn[-1] = 1800
print(xn)
# xn = np.hstack((xn_0, xn_1, xn_2, xn_3))
# # xn = xm + 300 #np.arange(-950.0, 1100.0, dx_rx)
yn = 0
zn = 1000
X = np.arange(-1400.0, 1800.0, dx)
Y = np.ones_like(X) - 1 # np.arange(-1000.0, 1000.0, dx)
Z = np.arange(0.0, 1025.0, dx)
x, z = np.meshgrid(X, Z)
y = 0.0
J1 = np.zeros((82, 256))
for idx in range(len(xm)):
  if xm[idx] != xa and xn[idx] != xa:
    # compactify the calculation
    arg1 = (((x - xa) * (x - xm[idx]) + (y - ya) * (y - ym) +
            (z - za) * (z - zm)) *
            (((x - xa)**2 + (y - ya)**2 + (z - za)**2)**-1.5) *
            (((x - xm[idx])**2 + (y - ym)**2 + (z - zm)**2)**-1.5))
    arg2 = (((x - xa) * (x - xn[idx]) + (y - ya) * (y - yn) +
            (z - za) * (z - zn)) *
            (((x - xa)**2 + (y - ya)**2 + (z - za)**2)**-1.5) *
            (((x - xn[idx])**2 + (y - yn)**2 + (z - zn)**2)**-1.5))
    arg3 = (((x - xb) * (x - xm[idx]) + (y - yb) * (y - ym) +
            (z - zb) * (z - zm)) *
            (((x - xb)**2 + (y - yb)**2 + (z - zb)**2)**-1.5) *
            (((x - xm[idx])**2 + (y - ym)**2 + (z - zm)**2)**-1.5))
    arg4 = (((x - xb) * (x - xn[idx]) + (y - yb) * (y - yn) +
            (z - zb) * (z - zn)) *
            (((x - xb)**2 + (y - yb)**2 + (z - zb)**2)**-1.5) *
            (((x - xn[idx])**2 + (y - yn)**2 + (z - zn)**2)**-1.5))
    # evaluate sensitivity equation
    J1 += (1 / (4 * np.pi**2) * (arg1 - arg2 - arg3 + arg4))**2
J1 = np.sqrt(J1 / len(xm))
ym_ = np.ones_like(xm) + zm
yn_ = np.ones_like(xn) + zn
ax = plt.subplot(5, 1, 1, aspect='equal')
ax.contourf(x, z, np.log(J1), cmap='viridis', vmin=-37, vmax=-25)
# pcolorOpts = {}
# ax.pcolormesh(x, z, np.log(J1), cmap='viridis', clim=(-32, -20), vmin=-32, vmax=-20, **pcolorOpts)
# pcolorOpts = {}
# ph = ax.pcolormesh(x, z, np.log(J1), cmap='viridis', clim=(-14, -32), vmin=-14, vmax=-32, **pcolorOpts)
# 
# ax.colorbar()
tx1 = [xa, za]
ax.plot(xm, ym_, 'k*')
ax.plot(xn, yn_, 'k*')
ax.plot(tx1[0], tx1[1], 'r*')
ax.grid()
ax.axes.set_xlim(-1400, 1800)
ax.set_title("1) expanding dipoles (a=50,100,200,400), 2) a=50, 3) a=300 4) a=300 rx-sep=50 5) Comp. Dipoles")
# plt.show()

# ============================================================================
# a=50
dx_rx = 50.0
xb = -6000.0
yb = 0.0
zb = 1000.0
xa = -1400.0
ya = 0.0
za = 1000.0
xm = np.arange(-1400.0, 1850.0, dx_rx)
ym = 0.0
zm = 1000.0
xn = xm + dx_rx #np.arange(-950.0, 1100.0, dx_rx)
yn = 0
zn = 1000
X = np.arange(-1400.0, 1800.0, dx)
Y = np.ones_like(X) - 1 # np.arange(-1000.0, 1000.0, dx)
Z = np.arange(0.0, 1025.0, dx)
x, z = np.meshgrid(X, Z)
y = 0.0
J2 = np.zeros((82, 256))
for idx in range(len(xm)):
  if xm[idx] != xa and xn[idx] != xa:
    # compactify the calculation
    arg1 = (((x - xa) * (x - xm[idx]) + (y - ya) * (y - ym) +
            (z - za) * (z - zm)) *
            (((x - xa)**2 + (y - ya)**2 + (z - za)**2)**-1.5) *
            (((x - xm[idx])**2 + (y - ym)**2 + (z - zm)**2)**-1.5))
    arg2 = (((x - xa) * (x - xn[idx]) + (y - ya) * (y - yn) +
            (z - za) * (z - zn)) *
            (((x - xa)**2 + (y - ya)**2 + (z - za)**2)**-1.5) *
            (((x - xn[idx])**2 + (y - yn)**2 + (z - zn)**2)**-1.5))
    arg3 = (((x - xb) * (x - xm[idx]) + (y - yb) * (y - ym) +
            (z - zb) * (z - zm)) *
            (((x - xb)**2 + (y - yb)**2 + (z - zb)**2)**-1.5) *
            (((x - xm[idx])**2 + (y - ym)**2 + (z - zm)**2)**-1.5))
    arg4 = (((x - xb) * (x - xn[idx]) + (y - yb) * (y - yn) +
            (z - zb) * (z - zn)) *
            (((x - xb)**2 + (y - yb)**2 + (z - zb)**2)**-1.5) *
            (((x - xn[idx])**2 + (y - yn)**2 + (z - zn)**2)**-1.5))
    # evaluate sensitivity equation
    J2 += (1 / (4 * np.pi**2) * (arg1 - arg2 - arg3 + arg4))**2
J2 = np.sqrt(J2 / len(xm))
ym_ = np.ones_like(xm) + zm
yn_ = np.ones_like(xn) + zn
ax2 = plt.subplot(5, 1, 2, aspect='equal')
ax2.contourf(x, z, np.log(J2), cmap='viridis', vmin=-37, vmax=-25)
tx1 = [xa, za]
ax2.plot(xm, ym_, 'k*')
ax2.plot(xn, yn_, 'k*')
ax2.plot(tx1[0], tx1[1], 'r*')
ax2.grid()
ax2.axes.set_xlim(-1400, 1800)
# ax2.set_title("a=50 Forward Dipoles")

# ============================================================================
# a=300
dx_rx = 300.0
xb = -6000.0
yb = 0.0
zb = 1000.0
xa = -1400.0
ya = 0.0
za = 1000.0
xm = np.arange(-1400.0, 1800, dx_rx) + 100
ym = 0.0
zm = 1000.0
xn = xm + dx_rx  #np.arange(-950.0, 1100.0, dx_rx)
yn = 0
zn = 1000
X = np.arange(-1400.0, 1800.0, dx)
Y = np.ones_like(X) - 1  # np.arange(-1000.0, 1000.0, dx)
Z = np.arange(0.0, 1025.0, dx)
x, z = np.meshgrid(X, Z)
y = 0.0
J = np.zeros((82, 256))
for idx in range(len(xm)):
  if xm[idx] != xa and xn[idx] != xa:
    # compactify the calculation
    arg1 = (((x - xa) * (x - xm[idx]) + (y - ya) * (y - ym) +
            (z - za) * (z - zm)) *
            (((x - xa)**2 + (y - ya)**2 + (z - za)**2)**-1.5) *
            (((x - xm[idx])**2 + (y - ym)**2 + (z - zm)**2)**-1.5))
    arg2 = (((x - xa) * (x - xn[idx]) + (y - ya) * (y - yn) +
            (z - za) * (z - zn)) *
            (((x - xa)**2 + (y - ya)**2 + (z - za)**2)**-1.5) *
            (((x - xn[idx])**2 + (y - yn)**2 + (z - zn)**2)**-1.5))
    arg3 = (((x - xb) * (x - xm[idx]) + (y - yb) * (y - ym) +
            (z - zb) * (z - zm)) *
            (((x - xb)**2 + (y - yb)**2 + (z - zb)**2)**-1.5) *
            (((x - xm[idx])**2 + (y - ym)**2 + (z - zm)**2)**-1.5))
    arg4 = (((x - xb) * (x - xn[idx]) + (y - yb) * (y - yn) +
            (z - zb) * (z - zn)) *
            (((x - xb)**2 + (y - yb)**2 + (z - zb)**2)**-1.5) *
            (((x - xn[idx])**2 + (y - yn)**2 + (z - zn)**2)**-1.5))
    # evaluate sensitivity equation
    J += (1 / (4 * np.pi**2) * (arg1 - arg2 - arg3 + arg4))**2
J = np.sqrt(J / len(xm))
ym_ = np.ones_like(xm) + zm
yn_ = np.ones_like(xn) + zn
ax3 = plt.subplot(5, 1, 3, aspect='equal')
ax3.contourf(x, z, np.log(J), cmap='viridis', vmin=-37, vmax=-25)
tx1 = [xa, za]
ax3.plot(xm, ym_, 'k*')
ax3.plot(xn, yn_, 'k*')
ax3.plot(tx1[0], tx1[1], 'r*')
ax3.grid()
ax3.axes.set_xlim(-1400, 1800)
# ax3.set_title("a=300 Forward and Reverse Dipoles")
# plt.show()

# ============================================================================
# a=300 rx-sep = 50

dx_rx = 50.0
xb = -6000.0
yb = 0.0
zb = 1000.0
xa = -1400.0
ya = 0.0
za = 1000.0
xm = np.arange(-1400.0, 1800, dx_rx)
ym = 0.0
zm = 1000.0
xn = xm + 300 #np.arange(-950.0, 1100.0, dx_rx)
print(xm, xn)
yn = 0
zn = 1000
X = np.arange(-1400.0, 1800.0, dx)
Y = np.ones_like(X) - 1 # np.arange(-1000.0, 1000.0, dx)
Z = np.arange(0.0, 1025.0, dx)
x, z = np.meshgrid(X, Z)
y = 0.0
J = np.zeros((82, 256))
for idx in range(len(xm)):
  if xm[idx] != xa and xn[idx] != xa:
    # compactify the calculation
    arg1 = (((x - xa) * (x - xm[idx]) + (y - ya) * (y - ym) +
            (z - za) * (z - zm)) *
            (((x - xa)**2 + (y - ya)**2 + (z - za)**2)**-1.5) *
            (((x - xm[idx])**2 + (y - ym)**2 + (z - zm)**2)**-1.5))
    arg2 = (((x - xa) * (x - xn[idx]) + (y - ya) * (y - yn) +
            (z - za) * (z - zn)) *
            (((x - xa)**2 + (y - ya)**2 + (z - za)**2)**-1.5) *
            (((x - xn[idx])**2 + (y - yn)**2 + (z - zn)**2)**-1.5))
    arg3 = (((x - xb) * (x - xm[idx]) + (y - yb) * (y - ym) +
            (z - zb) * (z - zm)) *
            (((x - xb)**2 + (y - yb)**2 + (z - zb)**2)**-1.5) *
            (((x - xm[idx])**2 + (y - ym)**2 + (z - zm)**2)**-1.5))
    arg4 = (((x - xb) * (x - xn[idx]) + (y - yb) * (y - yn) +
            (z - zb) * (z - zn)) *
            (((x - xb)**2 + (y - yb)**2 + (z - zb)**2)**-1.5) *
            (((x - xn[idx])**2 + (y - yn)**2 + (z - zn)**2)**-1.5))
    # evaluate sensitivity equation
    J += (1 / (4 * np.pi**2) * (arg1 - arg2 - arg3 + arg4))**2
J = np.sqrt(J / len(xm))
ym_ = np.ones_like(xm) + zm
yn_ = np.ones_like(xn) + zn
# create an axes to plot on
ax4 = plt.subplot(5, 1, 4, aspect='equal')
ax4.contourf(x, z, np.log(J), cmap='viridis', vmin=-37, vmax=-25)
tx1 = [xa, za]
ax4.plot(xm, ym_, 'k*')
ax4.plot(xn, yn_, 'k*')
ax4.plot(tx1[0], tx1[1], 'r*')
ax4.grid()
ax4.axes.set_xlim(-1400, 1800)
# ax4.set_title("a=300 rx-sep=50 Forward and Reverse Dipoles")
# plt.show()

# ============================================================================
# comp

dx_rx0 = 50.0
dx_rx1 = 100.0
dx_rx2 = 150.0
dx_rx3 = 200.0
dx_rx4 = 250.0
dx_rx5 = 300.0
dx_rx6 = 350.0
dx_rx7 = 400.0
dx_rx8 = 450.0
dx_rx9 = 550.0
dx_rx10 = 600.0
dx_rx11 = 650.0
dx_rx12 = 700.0
dx_rx13 = 750.0
dx_rx14 = 800.0
xb = -6000.0
yb = 0.0
zb = 1000.0
xa = -1400.0
ya = 0.0
za = 1000.0
xm0 = np.arange(-1400.0, 1800, dx_rx0)
xm1 = np.arange(-1400.0, 1800, dx_rx1)
xm2 = np.arange(-1400.0, 1800, dx_rx2)
xm3 = np.arange(-1400.0, 1800, dx_rx3)
xm4 = np.arange(-1400.0, 1800, dx_rx4)
xm5 = np.arange(-1400.0, 1800, dx_rx5)
xm6 = np.arange(-1400.0, 1800, dx_rx6)
xm7 = np.arange(-1400.0, 1800, dx_rx7)
xm8 = np.arange(-1400.0, 1800, dx_rx8)
xm9 = np.arange(-1400.0, 1800, dx_rx9)
xm10 = np.arange(-1400.0, 1800, dx_rx10)
xm11 = np.arange(-1400.0, 1800, dx_rx11)
xm12 = np.arange(-1400.0, 1800, dx_rx12)
xm13 = np.arange(-1400.0, 1800, dx_rx13)
xm14 = np.arange(-1400.0, 1800, dx_rx14)
xm = [xm0, xm1, xm2,xm3,xm4,xm5,xm6,xm7,xm8,xm9,xm10,xm11,xm12,xm13,xm14]
ym = 0.0
zm = 1000.0
xn = [xm0 + dx_rx0,xm1 + dx_rx1,xm2 + dx_rx2,xm3 + dx_rx3,xm4 + dx_rx4,xm5 + dx_rx5,xm6 + dx_rx6,xm7 + dx_rx7,xm8 + dx_rx8,xm9 + dx_rx9,xm10 + dx_rx10,xm11 + dx_rx11,xm12 + dx_rx12,xm13 + dx_rx13,xm14 + dx_rx14] 
# print(xm, xn)
yn = 0
zn = 1000
X = np.arange(-1400.0, 1800.0, dx)
Y = np.ones_like(X) - 1 # np.arange(-1000.0, 1000.0, dx)
Z = np.arange(0.0, 1025.0, dx)
x, z = np.meshgrid(X, Z)
y = 0.0
J = np.zeros((82, 256))
for idz in range(len(xm)):
  for idx in range(len(xm[idz])):
    if xm[idz][idx] != xa and xn[idz][idx] != xa:
      # compactify the calculation
      arg1 = (((x - xa) * (x - xm[idz][idx]) + (y - ya) * (y - ym) +
              (z - za) * (z - zm)) *
              (((x - xa)**2 + (y - ya)**2 + (z - za)**2)**-1.5) *
              (((x - xm[idz][idx])**2 + (y - ym)**2 + (z - zm)**2)**-1.5))
      arg2 = (((x - xa) * (x - xn[idz][idx]) + (y - ya) * (y - yn) +
              (z - za) * (z - zn)) *
              (((x - xa)**2 + (y - ya)**2 + (z - za)**2)**-1.5) *
              (((x - xn[idz][idx])**2 + (y - yn)**2 + (z - zn)**2)**-1.5))
      arg3 = (((x - xb) * (x - xm[idz][idx]) + (y - yb) * (y - ym) +
              (z - zb) * (z - zm)) *
              (((x - xb)**2 + (y - yb)**2 + (z - zb)**2)**-1.5) *
              (((x - xm[idz][idx])**2 + (y - ym)**2 + (z - zm)**2)**-1.5))
      arg4 = (((x - xb) * (x - xn[idz][idx]) + (y - yb) * (y - yn) +
              (z - zb) * (z - zn)) *
              (((x - xb)**2 + (y - yb)**2 + (z - zb)**2)**-1.5) *
              (((x - xn[idz][idx])**2 + (y - yn)**2 + (z - zn)**2)**-1.5))
      # evaluate sensitivity equation
      J += (1 / (4 * np.pi**2) * (arg1 - arg2 - arg3 + arg4))**2
J = np.sqrt(J / len(xm))
ax5 = plt.subplot(5, 1, 5, aspect='equal')
ax5.contourf(x, z, np.log(J), cmap='viridis', vmin=-37, vmax=-25)
tx1 = [xa, za]
# ax4.plot(xm, ym_, 'k*')
# ax4.plot(xn, yn_, 'k*')
ax5.plot(tx1[0], tx1[1], 'r*')
ax5.grid()
ax5.axes.set_xlim(-1400, 1800)
# ax5.set_title("comprehensive Forward and Reverse Dipoles")
plt.show()

# diff ===========================================================
plt.contourf(x, z, np.log(np.abs(J1)), cmap='viridis')
plt.colorbar()
plt.axes().set_aspect('equal', 'datalim')
plt.show()
