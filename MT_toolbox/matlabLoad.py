import matplotlib.pyplot as plt
import numpy as np
import DCIPtools as DCIP
import h5py
"""
   NOTES
   Matlab database contains line variables which has
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
output_path = 'E:/Projects/Qmag/IPHT'
lines_to_extract = 2

# load the database
f = h5py.File(file_path)

# USE THIS TO GET THE KEYS AND SUB ATTRIBUTES
# print(list(f.keys()))  # this what gave me the Line30 path
# f['Line30'].visititems(lambda n,o:print(n, o))

# get the data id
lines = list(f.keys())[1]

# =============================================================
# Lat long data retrieve ====================================

# Id tag for latitude
lati_id = lines + '/lat'
# Id tag for longitude
longi_id = lines + '/lat'
# altitude id tage
alti_id = lines + '/alt'

# get the latitude data
lat_convert = f.get(lati_id)
# get longitude
lon_convert = f.get(longi_id)
# get longitude
alt_convert = f.get(alti_id)
# create list to hold the data
lon_list = []
lat_list = []
alt_list = []
# retrieve data
for idx in range(lines_to_extract):
    lon_id = lon_convert[idx][0]
    lon_list.append(f[lon_id][0])
    lat_id = lat_convert[idx][0]
    lat_list.append(f[lat_id][0])
    alt_id = alt_convert[idx][0]
    alt_list.append(f[alt_id][0])


# ==========================================================
# get the reduced rotated data =============================

# mag data Id
mag_red_rot_id = lines + '/mag_temp_red_rot'
mag_reduced_rotated_convert = f.get(mag_red_rot_id)
# # create a nD X 3
mag_reduced_rotated = []
for idx in range(lines_to_extract):
    mag_reduced_rotated_id = mag_reduced_rotated_convert[idx][0]
    mag_reduced_rotated.append(f[mag_reduced_rotated_id])

# get high gain mag
# mag_high_gain_convert = f.get('Line30/mag_high')
# mag_high_gain_id = mag_high_gain_convert[0][0]
# # create a nD X 3
# mag_high_gain = f[mag_high_gain_id]

# # get low gain mag
# mag_low_gain_convert = f.get('Line30/mag_low')
# mag_low_gain_id = mag_low_gain_convert[0][0]
# # create a nD X 3
# mag_low_gain = f[mag_low_gain_id]

# # get mag temp?
# mag_temp_convert = f.get('Line30/mag_temp')
# mag_temp_id = mag_temp_convert[0][0]
# # create a nD X 3
# mag_temp = f[mag_temp_id]
# # create a nD X 3
# mag_temp = f[mag_temp_id]

# =========================================================
# retrieve IMU data =======================================

# gyros DAS id
gyros_das_id = lines + '/gyros_DAS'
# accelerations DAS id
acc_das_id = lines + '/acc_DAS'
# DAS attitudes id
att_das_id = lines + '/att_DAS'
# gyros uIMU id
gyros_imu_id = lines + '/gyros_uimu'
# accelerations uIMU id
acc_imu_id = lines + '/acc_uimu'
# gyros DAS id
att_imu_id = lines + '/att_DAS'

gyros_das_convert = f.get(gyros_das_id)
acc_das_convert = f.get(acc_das_id)
att_das_convert = f.get(att_das_id)
gyros_imu_convert = f.get(gyros_imu_id)
acc_imu_convert = f.get(acc_imu_id)
att_imu_convert = f.get(att_imu_id)

gyro_das = []
acc_das = []
att_das = []
gyro_imu = []
acc_imu = []
att_imu = []

for idx in range(lines_to_extract):
    gyro_das_id = gyros_das_convert[idx][0]
    gyro_das.append(f[gyro_das_id])
    acceleratrions_das_id = acc_das_convert[idx][0]
    acc_das.append(f[acceleratrions_das_id])
    attitude_das_id = att_das_convert[idx][0]
    att_das.append(f[attitude_das_id])

    gyro_imu_id = gyros_imu_convert[idx][0]
    gyro_imu.append(f[gyro_das_id])
    acceleratrions_imu_id = acc_imu_convert[idx][0]
    acc_imu.append(f[acceleratrions_imu_id])
    attitude_imu_id = att_imu_convert[idx][0]
    att_imu.append(f[attitude_imu_id])

# ======================================================================
# write a files =========================================================
for idx in range(lines_to_extract):
    lat_out = np.asarray(lat_list[idx])
    lat_path = output_path + '/lat' + str(idx + 1) + '.npy'
    np.save(lat_path, lat_out)
    lon_path = output_path + '/lon' + str(idx + 1) + '.npy'
    lon_out = np.asarray(lon_list[idx])
    np.save(lon_path, lon_out)
    alt_path = output_path + '/alt' + str(idx + 1) + '.npy'
    alt_out = np.asarray(alt_list[idx])
    np.save(alt_path, alt_out)
    data_path = output_path + '/mag_data_line' + str(idx + 1) + '.npy'
    data_out = np.asarray(mag_reduced_rotated[idx])
    np.save(data_path, data_out)

    gyro_das_path = output_path + '/gyro_das_line' + str(idx + 1) + '.npy'
    gyro_das_out = np.asarray(gyro_das[idx])
    np.save(gyro_das_path, gyro_das_out)
    acc_das_path = output_path + '/acc_das_line' + str(idx + 1) + '.npy'
    acc_das_out = np.asarray(acc_das[idx])
    np.save(acc_das_path, acc_das_out)
    att_das_path = output_path + '/att_das_line' + str(idx + 1) + '.npy'
    att_das_out = np.asarray(att_das[idx])
    np.save(att_das_path, att_das_out)

    gyro_imu_path = output_path + '/gyro_imu_line' + str(idx + 1) + '.npy'
    gyro_imu_out = np.asarray(gyro_imu[idx])
    np.save(gyro_imu_path, gyro_imu_out)
    acc_imu_path = output_path + '/acc_imu_line' + str(idx + 1) + '.npy'
    acc_imu_out = np.asarray(acc_imu[idx])
    np.save(acc_imu_path, acc_imu_out)
    att_imu_path = output_path + '/att_imu_line' + str(idx + 1) + '.npy'
    att_imu_out = np.asarray(att_imu[idx])
    np.save(att_imu_path, att_imu_out)

# Plotting ========================================
# print(len(acc_imu[0][0]))
# for idx in range(lines_to_extract):
#     # plt.plot(lon_list[idx], lat_list[idx])
#     plt.plot(acc_imu[0][0], acc_imu[0][1], 'ro')
# plt.title("coordinates (long vs lat)")
# plt.show()
