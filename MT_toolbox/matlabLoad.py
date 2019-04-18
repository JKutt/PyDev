import matplotlib.pyplot as plt
import numpy as np
import DCIPtools as DCIP
import h5py
"""
   NOTES
   Matlab database contains variables Line30 which has
   following attributes:

   acc_DAS
   acc_uimu
   alt
   att_DAS
   att_uimu
   lat
   lon
   mag_high
   mag_low
   mag_temp
   mag_temp_red_rot
"""

# filepath to the matlab database
file_path = "E:/Projects/Qmag/IPHT/Line30_20190328_211621.mat"

# load the database
f = h5py.File(file_path)

# USE THIS TO GET THE KEYS AND SUB ATTRIBUTES
# print(list(f.keys()))  # this what gave me the Line30 path
# f['Line30'].visititems(lambda n,o:print(n, o))

# =============================================================
# weird loading data stuff =====================================

# get the latitude data
lat_convert = f.get('Line30/lat')
lat_id = lat_convert[0][0]
lat = f[lat_id][0]
# get longitude data
lon_convert = f.get('Line30/lon')
lon_id = lon_convert[0][0]
lon = f[lon_id][0]

# get the reduced rotated data
mag_reduced_rotated_convert = f.get('Line30/mag_temp_red_rot')
mag_reduced_rotated_id = mag_reduced_rotated_convert[0][0]
# create a nD X 3
mag_reduced_rotated = f[mag_reduced_rotated_id]

# get high gain mag
mag_high_gain_convert = f.get('Line30/mag_high')
mag_high_gain_id = mag_high_gain_convert[0][0]
# create a nD X 3
mag_high_gain = f[mag_high_gain_id]

# get low gain mag
mag_low_gain_convert = f.get('Line30/mag_low')
mag_low_gain_id = mag_low_gain_convert[0][0]
# create a nD X 3
mag_low_gain = f[mag_low_gain_id]

# get mag temp?
mag_temp_convert = f.get('Line30/mag_temp')
mag_temp_id = mag_temp_convert[0][0]
# create a nD X 3
mag_temp = f[mag_temp_id]

plt.plot(lon, lat)
plt.title("coordinates (long vs lat)")
plt.show()

fig = plt.figure()
ax1 = fig.add_subplot(311)
ax2 = fig.add_subplot(312, sharex=ax1, sharey=ax1)
ax3 = fig.add_subplot(313, sharex=ax1, sharey=ax1)
ax1.set_title("X")
ax2.set_title("Y")
ax3.set_title("Z-")

ax1.plot(mag_reduced_rotated[0], 'r')
ax2.plot(mag_reduced_rotated[1], 'g')
ax3.plot(mag_reduced_rotated[2], 'b')

# ax1.plot(mag_high_gain[0], 'k')
# ax2.plot(mag_high_gain[1], 'c')
# ax3.plot(mag_high_gain[2], 'm')

# ax1.plot(mag_low_gain[0], '-.k')
# ax2.plot(mag_low_gain[1], '-.c')
# ax3.plot(mag_low_gain[2], '-.m')

# ax1.plot(mag_temp[0], '.k')
# ax2.plot(mag_temp[1], '.c')
# ax3.plot(mag_temp[2], '.m')

plt.show()

# padd = DCIP.padNextPower2(np.asarray(mag_reduced_rotated[1]))
# padd_fft = DCIP.getSpectralDensity(padd)
# freqs = np.arange(0, padd_fft.size) * (31250. / padd_fft.size)
# plt.semilogy(freqs[0:int(padd_fft.size/2)], padd_fft[0:int(padd_fft.size/2)])
# plt.show()
# write a file
# lat_out = np.asarray(lat)
# np.save('E:/Projects/Qmag/IPHT/lat.npy', lat_out)
# lon_out = np.asarray(lon)
# np.save('E:/Projects/Qmag/IPHT/lon.npy', lon_out)