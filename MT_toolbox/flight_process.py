import numpy as np
import utm
import pylab as plt
import scipy.signal as scygnal
import glob


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
    fOut.write('<processing_parameters>\n')
    fOut.write('\t<tbw=' + args[12] + '>\n')
    fOut.write('\t<nfft=' + args[13] + '>\n')
    fOut.write('\t<overlap=' + args[14] + '>\n')
    fOut.write('\t<taper=' + args[15] + '>\n')
    fOut.write('\t<index_first_frequency=' + args[16] + '>\n')
    fOut.write('\t<frequency_increment=' + args[17] + '>\n')
    fOut.write('\t<nb_increment=' + args[18] + '>\n')
    fOut.write('\t<length_reduction=' + args[19] + '>\n')
    fOut.write('\t<nb_reductions=' + args[20] + '>\n')
    fOut.write('\t<output=' + args[21] + '>\n')
    fOut.write('\t<error=' + args[22] + '>\n')
    fOut.write('\t<remote_coherence=' + args[23] + '>\n')

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

path = 'C:/Users/HugoLarnier/Desktop/Data/flight_test/20190328/'

lat = np.load(path + 'lat.npy')
lon = np.load(path + 'lon.npy')

if not len(glob.glob('utm_flight.npy')):
    utm_flight = np.zeros((len(lat), 2))
    for latitude, longitude, i in zip(lat, lon, range(len(lat))):
        print('Percentage: ', str(i / len(lat) * 100.), '%')
        tmp = utm.from_latlon(latitude, longitude)
        utm_flight[i, 0] = tmp[0]
        utm_flight[i, 1] = tmp[1]

    utm_flight = np.asarray(utm_flight)
    np.save('utm_flight.npy', utm_flight)
else:
    utm_flight = np.load('utm_flight.npy')

data = np.load(path + 'mag_rotated_reduced.npy')

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


nbStations = 20
nfft = 131072
for i in range(nbStations):
    fOut = open('./processing/hx' + str(i + 1) + '.dat', 'w')
    for j in range(nfft):
        fOut.write(str(X[i * nfft + j]) + '\n')
    fOut.close()

    fOut = open('./processing/hy' + str(i + 1) + '.dat', 'w')
    for j in range(131072):
        fOut.write(str(Y[i * nfft + j]) + '\n')
    fOut.close()

    fOut = open('./processing/hz' + str(i + 1) + '.dat', 'w')
    for j in range(131072):
        fOut.write(str(Z[i * nfft + j]) + '\n')
    fOut.close()
    path_proc = 'C:/Users/HugoLarnier/Desktop/Projects/QMAG/flight_20190328/dev/processing/'
    arguments = ['Flight test 20190328',
                 'MT_' + str(i + 1),
                 'DIAS',
                 'Hugo Larnier',
                 'C:/Users/HugoLarnier/Desktop/Projects/QMAG/flight_20190328/dev/processing/',
                 'T',
                 '1',
                 path_proc + 'hz' + str(i + 1) + '.dat',
                 path_proc + 'hx' + str(i + 1) + '.dat',
                 path_proc + 'hy' + str(i + 1) + '.dat',
                 path_proc + 'hx' + str(i + 1) + '.dat',
                 path_proc + 'hy' + str(i + 1) + '.dat',
                 '2',
                 '1024',
                 '20',
                 'slepian',
                 '8',
                 '2',
                 '4',
                 '2',
                 '5',
                 'MT_' + str(i + 1),
                 'jackknife',
                 '0.0']
    write_inp_file('C:/Users/HugoLarnier/Desktop/Projects/QMAG/flight_20190328/dev/processing/', arguments)
    # Call to diasMT here


exit()
fig, ax = plt.subplots(2, 1, sharex=True)
ax[0].plot(utm_flight[:, 1])
ax[1].plot(data[0])
plt.show()
fig, ax = plt.subplots(3, 1, sharex=True)
for axes, i in zip(ax, range(len(ax))):
    axes.plot(time, data[i])
plt.show()

fig, ax = plt.subplots(3, 1, sharex=True)
for axes, i in zip(ax, range(len(ax))):
    axes.set_xscale('log')
    axes.set_yscale('log')
    axes.plot(freq, np.abs(np.fft.fft(data[i])))
plt.show()
