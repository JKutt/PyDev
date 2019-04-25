import numpy as np
import utm
import pylab as plt
import scipy.signal as scygnal
import glob
import matplotlib.cm as cm
from scipy.interpolate import griddata

plt.rcParams.update({'font.size': 22})
def read_mt_station(path):

    nHeader = 14
    fIn = open(path, 'r')
    lines = fIn.readlines()
    nFreq = len(lines) - nHeader
    frequencies = [0. for i in range(nFreq)]
    for (line, i) in zip(lines, range(len(lines))):
        if '#Location' in line:
            tmp = line.split('#Location (UTM):')[1].split(';')
            location = [float(tmp[0]), float(tmp[1]), float(tmp[2])]

        if '#Tensor type:' in line:
            type = line.split('#Tensor type:')[1][1:-1]
            if type == 'Tipper':
                nChannels = 1
            if type == 'Z':
                nChannels = 2
            if type == 'Full':
                nChannels = 3

            Z = np.asarray([[[np.complex(0, 0), np.complex(0, 0)]
                                         for f in range(nFreq)]
                                         for channel in range(nChannels)])
            error_Z = np.asarray([[[0., 0.]
                                         for f in range(nFreq)]
                                         for channel in range(nChannels)])

        if i >= nHeader:
            frequencies[i - nHeader] = float(line[0:15])
            for channel in range(nChannels):
                for j in range(2):
                    Z[channel, i - nHeader, j] += float(line[(j + 1 + channel) * 15: (j + 2 + channel) * 15])
            end = (j + 2 + channel) * 15
            for channel in range(nChannels):
                for j in range(2):
                    Z[channel, i - nHeader, j] += float(line[(j + channel) * 15 + end:
                                                              (j + channel + 1) * 15 + end]) * np.complex(0, 1)
            end = (j + channel + 1) * 15 + end
            for channel in range(nChannels):
                for j in range(2):
                    error_Z[channel, i - nHeader, j] = float(line[(j + channel) * 15 + end:
                                                              (j + channel + 1) * 15 + end])

    return location, frequencies, Z, error_Z

path = 'C:/Users/HugoLarnier/Desktop/Projects/QMAG/flight_20190328/dev/processing/line1/nfft_65536/'
stations = [path + 'MT_' + str(i + 1) + '.DAT' for i in range(10)]

location, frequencies, Z, error_Z = read_mt_station(stations[0])
lines = ['line' + str(i + 1) for i in range(31)]
signY = [1, -1, 1, -1, 1, -1, 1, -1, 1 ,-1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1]
signX = [1 for i in range(len(lines))]
signY = signX[:]
#signX[18] *= -1
#signX[17] *= -1
#sign = [1 for i in range(len(stations))]
vmin = -0.5
vmax = 0.5

for ind_freq in range(2, len(frequencies)):

    X = []
    Y = []
    Tx = []
    Ty = []
    ETx = []
    ETy = []
    for line, ind in zip(lines, range(0, 13)):
        path = 'C:/Users/HugoLarnier/Desktop/Projects/QMAG/flight_20190328/dev/processing/' + line + '/'
        stations = glob.glob(path + 'MT*.DAT')
        for station in stations:
            print('Reading station: ', station)
            location, frequencies, Z, error_Z = read_mt_station(station)
            X.extend([location[1]])
            Y.extend([location[0]])
            Tx.extend([signX[ind] * Z[0, ind_freq, 0]])
            Ty.extend([signY[ind] * Z[0, ind_freq, 1]])
            #print(Ty[-1], Z[0, ind_freq, 1])
            ETx.extend([error_Z[0, ind_freq, 0]])
            ETy.extend([error_Z[0, ind_freq, 1]])

    Tx = np.asarray(Tx)
    Ty = np.asarray(Ty)
    X = np.asarray(X) / 1000
    Y = np.asarray(Y) / 1000

    fig, ((ax0, ax1, ax4), (ax2, ax3, ax5)) = plt.subplots(2, 3, sharex=True, figsize=(30, 20))
    ax0.set_xlabel('Easting')
    ax0.set_ylabel('Northing')
    ax1.set_xlabel('Easting')
    ax1.set_ylabel('Northing')
    fig.suptitle('Frequency:' + str(int(frequencies[ind_freq])) + ' Hz')

    Txclbr = ax0.scatter(Y, X, c=np.real(Tx), s=100, cmap="coolwarm", vmin=vmin, vmax=vmax)
    cbar = fig.colorbar(Txclbr, ax=ax0)
    cbar.ax.set_ylabel('Real Tx')

    Tyclbr = ax1.scatter(Y, X, c=np.imag(Tx), s=100, cmap="coolwarm", vmin=vmin, vmax=vmax)
    cbar = fig.colorbar(Tyclbr, ax=ax1)
    cbar.ax.set_ylabel('Imaginary Tx')

    Txclbr = ax2.scatter(Y, X, c=np.real(Ty), s=100, cmap="coolwarm", vmin=vmin, vmax=vmax)
    cbar = fig.colorbar(Txclbr, ax=ax2)
    cbar.ax.set_ylabel('Real Ty')

    Tyclbr = ax3.scatter(Y, X, c=np.imag(Ty), s=100, cmap="coolwarm", vmin=vmin, vmax=vmax)
    cbar = fig.colorbar(Tyclbr, ax=ax3)
    cbar.ax.set_ylabel('Imaginary Ty')


    Txclbr = ax4.scatter(Y, X, c=ETx, s=100, cmap="coolwarm", vmin=vmin, vmax=vmax)
    cbar = fig.colorbar(Txclbr, ax=ax4)
    cbar.ax.set_ylabel('Error Tx')

    Tyclbr = ax5.scatter(Y, X, c=ETy, s=100, cmap="coolwarm", vmin=vmin, vmax=vmax)
    cbar = fig.colorbar(Tyclbr, ax=ax5)
    cbar.ax.set_ylabel('Error Ty')

    X = []
    Y = []
    Tx = []
    Ty = []
    ETx = []
    ETy = []

    for line, ind in zip(lines, range(13, len(lines))):
        path = 'C:/Users/HugoLarnier/Desktop/Projects/QMAG/flight_20190328/dev/processing/' + lines[ind] + '/'
        stations = glob.glob(path + 'MT*.DAT')
        for station in stations:
            print('Reading station: ', station)
            location, frequencies, Z, error_Z = read_mt_station(station)
            X.extend([location[1]])
            Y.extend([location[0]])
            Tx.extend([signX[ind] * Z[0, ind_freq, 0]])
            Ty.extend([signY[ind] * Z[0, ind_freq, 1]])
            #print(Ty[-1], Z[0, ind_freq, 1])
            ETx.extend([error_Z[0, ind_freq, 0]])
            ETy.extend([error_Z[0, ind_freq, 1]])


    Tx = np.asarray(Tx)
    Ty = np.asarray(Ty)
    X = np.asarray(X) / 1000
    Y = np.asarray(Y) / 1000

    Txclbr = ax0.scatter(Y, X, c=np.real(Tx), s=100, cmap="coolwarm", vmin=vmin, vmax=vmax, marker='v')

    Tyclbr = ax1.scatter(Y, X, c=np.imag(Tx), s=100, cmap="coolwarm", vmin=vmin, vmax=vmax, marker='v')

    Txclbr = ax2.scatter(Y, X, c=np.real(Ty), s=100, cmap="coolwarm", vmin=vmin, vmax=vmax, marker='v')

    Tyclbr = ax3.scatter(Y, X, c=np.imag(Ty), s=100, cmap="coolwarm", vmin=vmin, vmax=vmax, marker='v')

    Txclbr = ax4.scatter(Y, X, c=ETx, s=100, cmap="coolwarm", vmin=vmin, vmax=vmax, marker='v')

    Tyclbr = ax5.scatter(Y, X, c=ETy, s=100, cmap="coolwarm", vmin=vmin, vmax=vmax, marker='v')

    Tx = np.asarray(Tx)
    Ty = np.asarray(Ty)
    X = np.asarray(X) / 1000
    Y = np.asarray(Y) / 1000

    plt.savefig('Map_' + str(ind_freq + 1) + '.png')
    #plt.show()

    # fig, ax0 = plt.subplots(1, 1, figsize=(30, 20))
    # for i in range(len(Y)):
    #     fig.suptitle('Frequency:' + str(int(frequencies[ind_freq])) + ' Hz')
    #
    #     scale = 0.05
    #     #ax0.arrow(10., 0., 0.5, 0.5, head_width=0.05, head_length=0.1, fc='k', ec='k')
    #     ax0.set_xlim([401, 402])
    #     ax0.set_ylim([5723, 5726])
    #     ax0.arrow(Y[i], X[i], -np.real(Tx[i]) * scale, -np.real(Ty[i]) * scale, width=0.001, head_width=0.01, head_length=0.01, fc='k', ec='k')
    plt.show()

exit()
for ind_freq in range(len(frequencies)):
    fig, (ax0, ax1) = plt.subplots(2, 1, sharex=True, figsize=(30, 20))
    ax0.set_xlabel('Northing')
    ax0.set_ylabel('Tx [ ]')
    ax1.set_xlabel('Northing')
    ax1.set_ylabel('Ty [ ]')

    X = []
    Tx = []
    Ty = []
    ETx = []
    ETy = []
    path = 'C:/Users/HugoLarnier/Desktop/Projects/QMAG/flight_20190328/dev/processing/line1/nfft_65536/'
    stations = [path + 'MT_' + str(i + 1) + '.DAT' for i in range(40)]
    for station in stations:
        print('Reading station: ', station)
        location, frequencies, Z, error_Z = read_mt_station(station)
        X.extend([location[1]])
        Tx.extend([Z[0, ind_freq, 0]])
        Ty.extend([Z[0, ind_freq, 1]])
        ETx.extend([error_Z[0, ind_freq, 0]])
        ETy.extend([error_Z[0, ind_freq, 1]])
    Tx = np.asarray(Tx)
    Ty = np.asarray(Ty)
    fig.suptitle('Frequency:' + str(int(frequencies[ind_freq])) + ' Hz')
    ax0.errorbar(X, np.real(Tx), yerr=ETx, fmt='rv', label='Real Tx, nfft=65536', markersize=12)
    ax0.errorbar(X, np.imag(Tx), yerr=ETx, fmt='bv', label='Imaginary Tx, nfft=65536', markersize=12)
    ax0.legend(loc="upper right")
    ax1.errorbar(X, np.real(Ty), yerr=ETy, fmt='rv', label='Real Ty, nfft=65536', markersize=12)
    ax1.errorbar(X, np.imag(Ty), yerr=ETy, fmt='bv', label='Imaginary Ty, nfft=65536', markersize=12)
    ax1.legend(loc="upper right")


    X = []
    Tx = []
    Ty = []
    ETx = []
    ETy = []
    path = 'C:/Users/HugoLarnier/Desktop/Projects/QMAG/flight_20190328/dev/processing/line2/nfft_65536/'
    stations = [path + 'MT_' + str(i + 1) + '.DAT' for i in range(40)]
    for station in stations:
        print('Reading station: ', station)
        location, frequencies, Z, error_Z = read_mt_station(station)
        X.extend([location[1]])
        Tx.extend([Z[0, ind_freq, 0]])
        Ty.extend([Z[0, ind_freq, 1]])
        ETx.extend([error_Z[0, ind_freq, 0]])
        ETy.extend([error_Z[0, ind_freq, 1]])
    Tx = np.asarray(Tx)
    Ty = np.asarray(Ty)
    fig.suptitle('Frequency:' + str(int(frequencies[ind_freq])) + ' Hz')
    ax0.errorbar(X, - np.real(Tx), yerr=ETx, fmt='ro', label='Real Tx, nfft=65536', markersize=12)
    ax0.errorbar(X, - np.imag(Tx), yerr=ETx, fmt='bo', label='Imaginary Tx, nfft=65536', markersize=12)
    ax0.legend(loc="upper right")
    ax1.errorbar(X, - np.real(Ty), yerr=ETy, fmt='ro', label='Real Ty, nfft=65536', markersize=12)
    ax1.errorbar(X, - np.imag(Ty), yerr=ETy, fmt='bo', label='Imaginary Ty, nfft=65536', markersize=12)
    ax1.legend(loc="upper right")

    X = []
    Tx = []
    Ty = []
    ETx = []
    ETy = []
    path = 'C:/Users/HugoLarnier/Desktop/Projects/QMAG/flight_20190328/dev/processing/line3/nfft_65536/'
    stations = [path + 'MT_' + str(i + 1) + '.DAT' for i in range(63)]
    for station in stations:
        print('Reading station: ', station)
        location, frequencies, Z, error_Z = read_mt_station(station)
        X.extend([location[1]])
        Tx.extend([Z[0, ind_freq, 0]])
        Ty.extend([Z[0, ind_freq, 1]])
        ETx.extend([error_Z[0, ind_freq, 0]])
        ETy.extend([error_Z[0, ind_freq, 1]])
    Tx = np.asarray(Tx)
    Ty = np.asarray(Ty)
    fig.suptitle('Frequency:' + str(int(frequencies[ind_freq])) + ' Hz')
    ax0.errorbar(X, np.real(Tx), yerr=ETx, fmt='rx', label='Real Tx, nfft=65536', markersize=12)
    ax0.errorbar(X, np.imag(Tx), yerr=ETx, fmt='bx', label='Imaginary Tx, nfft=65536', markersize=12)
    ax0.legend(loc="upper right")
    ax1.errorbar(X, np.real(Ty), yerr=ETy, fmt='rx', label='Real Ty, nfft=65536', markersize=12)
    ax1.errorbar(X, np.imag(Ty), yerr=ETy, fmt='bx', label='Imaginary Ty, nfft=65536', markersize=12)
    ax1.legend(loc="upper right")

    X = []
    Tx = []
    Ty = []
    ETx = []
    ETy = []
    path = 'C:/Users/HugoLarnier/Desktop/Projects/QMAG/flight_20190328/dev/processing/line4/nfft_65536/'
    stations = [path + 'MT_' + str(i + 1) + '.DAT' for i in range(42)]
    for station in stations:
        print('Reading station: ', station)
        location, frequencies, Z, error_Z = read_mt_station(station)
        X.extend([location[1]])
        Tx.extend([Z[0, ind_freq, 0]])
        Ty.extend([Z[0, ind_freq, 1]])
        ETx.extend([error_Z[0, ind_freq, 0]])
        ETy.extend([error_Z[0, ind_freq, 1]])
    Tx = np.asarray(Tx)
    Ty = np.asarray(Ty)
    fig.suptitle('Frequency:' + str(int(frequencies[ind_freq])) + ' Hz')
    ax0.errorbar(X, - np.real(Tx), yerr=ETx, fmt='rx', label='Real Tx, nfft=65536', markersize=12)
    ax0.errorbar(X, - np.imag(Tx), yerr=ETx, fmt='bx', label='Imaginary Tx, nfft=65536', markersize=12)
    ax0.legend(loc="upper right")
    ax1.errorbar(X, - np.real(Ty), yerr=ETy, fmt='rx', label='Real Ty, nfft=65536', markersize=12)
    ax1.errorbar(X, - np.imag(Ty), yerr=ETy, fmt='bx', label='Imaginary Ty, nfft=65536', markersize=12)
    ax1.legend(loc="upper right")
    plt.savefig('Profile_' + str(ind_freq + 1) + '.png')
    plt.show()
