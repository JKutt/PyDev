# Code for automatic processing of QMAG data from the '.mat' file
# Based on John Kuttai's loader

import matplotlib.pyplot as plt
import numpy as np
import math
import h5py
from glob import glob
from os import mkdir
import utm
import scipy.signal as scygnal
from subprocess import Popen, PIPE
from numba import vectorize, guvectorize

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

def write_inp_file(path, *args):

    args = args[0]
    fOut = open(path + '/parameters.inp', 'w')
    fOut.write('<project>\n')
    fOut.write('\t<patch=' + args[0] + '>\n')
    fOut.write('\t<station=' + args[1] + '>\n')
    fOut.write('\t<client=' + args[2] + '>\n')
    fOut.write('\t<processor=' + args[3] + '>\n')
    fOut.write('\t<folder=' + args[4] + '>\n')
    fOut.write('\t<type=' + args[5] + '>\n')
    fOut.write('\t<output_level=' + args[6] + '>\n')
    fOut.write('<data>\n')
    fOut.write('\t<hz=' + args[7] + '>\n')
    fOut.write('\t<hx=' + args[8] + '>\n')
    fOut.write('\t<hy=' + args[9] + '>\n')
    fOut.write('\t<rx=' + args[10] + '>\n')
    fOut.write('\t<ry=' + args[11] + '>\n')
    fOut.write('\t<location=' + args[12] + '>\n')
    fOut.write('<processing_parameters>\n')
    fOut.write('\t<tbw=' + args[13] + '>\n')
    fOut.write('\t<nfft=' + args[14] + '>\n')
    fOut.write('\t<overlap=' + args[15] + '>\n')
    fOut.write('\t<taper=' + args[16] + '>\n')
    fOut.write('\t<index_first_frequency=' + args[17] + '>\n')
    fOut.write('\t<frequency_increment=' + args[18] + '>\n')
    fOut.write('\t<nb_increment=' + args[19] + '>\n')
    fOut.write('\t<length_reduction=' + args[20] + '>\n')
    fOut.write('\t<nb_reductions=' + args[21] + '>\n')
    fOut.write('\t<output=' + args[22] + '>\n')
    fOut.write('\t<error=' + args[23] + '>\n')
    fOut.write('\t<remote_coherence=' + args[24] + '>\n')

    fOut.close()

## Function to filter low frequencies in VHF time series if they are too long
def HighPassFilter(data, sampleFreq, orderFilter, freqCut):
    nfft = len(data)
    b, a = scygnal.butter(orderFilter, np.asarray(freqCut / sampleFreq * 2.), 'high') #
    tmp = scygnal.filtfilt(b, a, np.asarray(data))
    data = tmp
    return data


def LowPassFilter(data, sampleFreq, orderFilter, freqCut):
    nfft = len(data)
    b, a = scygnal.butter(orderFilter, np.asarray(freqCut / sampleFreq * 2.), 'low') #
    tmp = scygnal.filtfilt(b, a, np.asarray(data))
    data = tmp
    return data


@guvectorize(['void(float32[:], float32[:])'], '(n)->(n)', target='cuda')
def deg2utm(p, q):
    size_of_data = int(p.size / 2)
    for idx in range(size_of_data):
        sa = 6378137.000000
        sb = 6356752.314245
        e2 = ((sa**2 - sb**2)**0.5) / sb
        e2cuadrada = e2**2
        c = sa**2 / sb
        lat = p[idx] * (np.pi / 180.0)
        lon = p[idx + size_of_data] * (np.pi / 180.0)
        Huso = math.floor((p[idx + size_of_data] / 6) + 31)
        S = ((Huso * 6) - 183)
        deltaS = lon - (S * (np.pi / 180))
        a = math.cos(lat) * math.sin(deltaS)
        epsilon = 0.5 * math.log((1 + a) / (1 - a))
        nu = math.atan(math.tan(lat) / math.cos(deltaS)) - lat
        v = (c / (1.0 + e2cuadrada * (math.cos(lat))**2)**0.5) * 0.9996
        ta = (e2cuadrada / 2.0) * epsilon**2 * (math.cos(lat))**2
        a1 = math.sin(2.0 * lat)
        a2 = a1 * (math.cos(lat))**2
        j2 = lat + (a1 / 2.0)
        j4 = ((3.0 * j2) + a2) / 4.0
        j6 = ((5.0 * j4) + (a2 * (math.cos(lat))**2)) / 3.0
        alfa = (3.0 / 4.0) * e2cuadrada
        beta = (5.0 / 3.0) * alfa**2
        gama = (35.0 / 27.0) * alfa**3
        Bm = 0.9996 * c * (lat - alfa * j2 + beta * j4 - gama * j6)
        q[idx] = epsilon * v * (1.0 + (ta / 3.0)) + 500000
        q[idx + size_of_data] = nu * v * (1.0 + ta) + Bm

        if q[idx + size_of_data] < 0:
            q[idx + size_of_data] = 9999999 + q[idx + size_of_data]


# Path to python
python_path = 'C:/ProgramData/Anaconda3/python.exe'
# Path to diasMT.py
diasMT_path = 'E:/afmag/PyDev/MT_toolbox/diasMT.py'
# filepath to the matlab database
file_path = "E:/afmag/data/Line30_20190328_211621.mat"
# folder path to extract .npy files
output_path = 'E:/afmag/flight_20190328/data'
# Time series length to create QMAG stations
nfft = int(131072 / 2)

lines_to_extract = 36

# load the database
f = h5py.File(file_path)

# USE THIS TO GET THE KEYS AND SUB ATTRIBUTES
#print(list(f.keys()))  # this what gave me the Line30 path
#f['Line30'].visititems(lambda n,o:print(n, o))
#exit()
# get the data id
lines = list(f.keys())[1]

lati_id = lines + '/lat'
# Id tag for longitude
longi_id = lines + '/lon'
# altitude id tage
alti_id = lines + '/alt'
# rotated mag id
mag_red_rot_id = lines + '/mag_temp_red_rot'
#mag_red_rot_id = lines + '/mag_temp'

# get the latitude data
lat_convert = f.get(lati_id)
# get longitude
lon_convert = f.get(longi_id)
# get longitude
alt_convert = f.get(alti_id)
# get mag
mag_reduced_rotated_convert = f.get(mag_red_rot_id)

for idx in range(lines_to_extract):
    print('Processing Line ' + str(idx))
    lon_id = lon_convert[idx][0]
    lon_list = f[lon_id][0]
    lat_id = lat_convert[idx][0]
    lat_list = f[lat_id][0]
    alt_id = alt_convert[idx][0]
    alt_list = f[alt_id][0]

    # # create a nD X 3

    mag_reduced_rotated_id = mag_reduced_rotated_convert[idx][0]
    mag_reduced_rotated = f[mag_reduced_rotated_id]

    tmp = glob(output_path + '/line' + str(idx))
    if not len(tmp):
        mkdir(output_path + '/line' + str(idx))
    else:
        print('Extraction folder already created.')

    lat_path = output_path + '/line' + str(idx) + '/lat.npy'
    np.save(lat_path, lat_list)
    lon_path = output_path + '/line' + str(idx) + '/lon.npy'
    np.save(lon_path, lon_list)
    alt_path = output_path + '/line' + str(idx) + '/alt.npy'
    np.save(alt_path, alt_list)
    data_path = output_path + '/line' + str(idx) + '/mag_data.npy'
    np.save(data_path, mag_reduced_rotated)

    lat = np.load(output_path + '/line' + str(idx) + '/lat.npy')
    lat = np.array(lat, dtype='f')
    lon = np.load(output_path + '/line' + str(idx) + '/lon.npy')
    lon = np.array(lon, dtype='f')
    safe_size = 100000
    chunks = int(np.floor(lat.size / safe_size))
    # print("data size: {0} & chunks {1}".format(lat.size, chunks))
    if not len(glob(output_path + '/line' + str(idx) + '/utm_flight.npy')):
        utm_flight = np.zeros((len(lat), 2))
        for idz in range(chunks + 1):
            # print("on chunk: ", idz)
            if idz == chunks:
                start = idz * safe_size
                end = lat.size
                v = np.hstack((lat[start:end],lon[start:end]))
                convert_out = deg2utm(v)
                size_ = int(convert_out.size / 2)
                utm_flight[start:end, 1] = convert_out[:size_]
                utm_flight[start:end, 0] = convert_out[size_:]
            else:
                start = idz * safe_size
                end = (idz + 1) * safe_size
                v = np.hstack((lat[start:end],lon[start:end]))
                convert_out = deg2utm(v)
                size_ = int(convert_out.size / 2)
                utm_flight[start:end, 1] = convert_out[:size_]
                utm_flight[start:end, 0] = convert_out[size_:]
        # for latitude, longitude, i in zip(lat, lon, range(len(lat))):
        #     print('UTM conversion percentage: ', str(i / len(lat) * 100.), '%')
        #     tmp = utm.from_latlon(latitude, longitude)
        #     utm_flight[i, 0] = tmp[0]
        #     utm_flight[i, 1] = tmp[1]

        utm_flight = np.asarray(utm_flight)
        np.save(output_path + '/line' + str(idx) + '/utm_flight.npy', utm_flight)
    else:
        utm_flight = np.load(output_path + '/line' + str(idx) + '/utm_flight.npy')

    data = np.load(output_path + '/line' + str(idx) + '/mag_data.npy')

    sample_freq = 31250
    X = data[0]
    Y = data[1]
    Z = data[2]

    #X = LowPassFilter(X, sample_freq, 1, 2000)
    #Y = LowPassFilter(Y, sample_freq, 1, 2000)
    #Z = LowPassFilter(Z, sample_freq, 1, 2000)

    X = HighPassFilter(X, sample_freq, 1, 10)
    Y = HighPassFilter(Y, sample_freq, 1, 10)
    Z = HighPassFilter(Z, sample_freq, 1, 10)

    data = [X, Y, Z]
    time = np.linspace(1, len(X), len(X)) / sample_freq
    freq = np.linspace(0, 1, len(X)) * sample_freq

    nbStations = len(data)
    nbStations = int(len(data[0]) / nfft) - 1

    for i in range(nbStations):
        path_proc = 'E:/afmag/flight_20190328/dev/processing/'
        tmp = glob(path_proc + '/line' + str(idx))
        if not len(tmp):
            mkdir(path_proc + '/line' + str(idx))
        else:
            print('Processing folder already created.')

        path_proc = 'E:/afmag/flight_20190328/dev/processing/line' + str(idx)
        fOut = open(path_proc + '/hx' + str(i + 1) + '.dat', 'w')
        for j in range(nfft):
            fOut.write(str(X[i * nfft + j]) + '\n')
        fOut.close()

        fOut = open(path_proc + '/hy' + str(i + 1) + '.dat', 'w')
        for j in range(nfft):
            fOut.write(str(Y[i * nfft + j]) + '\n')
        fOut.close()

        fOut = open(path_proc + '/hz' + str(i + 1) + '.dat', 'w')
        for j in range(nfft):
            fOut.write(str(Z[i * nfft + j]) + '\n')
        fOut.close()

        location = [float(np.median(utm_flight[i * nfft:(i + 1) * nfft, 0])),
                    float(np.median(utm_flight[i * nfft:(i + 1) * nfft, 1])),
                    0., '12']

        arguments = ['Flight test 20190328',
                     'MT_' + str(i + 1),
                     'DIAS',
                     'Hugo Larnier',
                     path_proc,
                     'T',
                     '1',
                     path_proc + '/hz' + str(i + 1) + '.dat',
                     path_proc + '/hx' + str(i + 1) + '.dat',
                     path_proc + '/hy' + str(i + 1) + '.dat',
                     path_proc + '/hx' + str(i + 1) + '.dat',
                     path_proc + '/hy' + str(i + 1) + '.dat',
                     str(location[0]) + ';' + str(location[1]) + ';' + str(location[2]) + ';Zone:' + location[3],
                     '2',
                     '1024',
                     '20',
                     'slepian',
                     '8',
                     '2',
                     '1',
                     '2',
                     '5',
                     'MT_' + str(i + 1),
                     'jackknife',
                     '0.0']
        write_inp_file('E:/afmag/PyDev/MT_toolbox/', arguments)
        # Call to diasMT here
        print('Starting MT process')
        print([python_path, diasMT_path])
        process = Popen([python_path, diasMT_path, '-c'], stdout=PIPE, stderr=PIPE, shell=True)

        #stdout = process.communicate()

        with process.stdout:
            for line in iter(process.stdout.readline, b''):
                print(line)

        process.wait()
        print('Ending MT process')

exit()
