#!C:\\Users\\HugoLarnier\\Anaconda3\\python

from datfiles_lib import *
import csv
import matplotlib.pylab as plt
import matplotlib.lines as mlines
from glob import glob
import numpy as np
from sys import exit
import matplotlib.cm as cm
import signalAnalysis as signal
import utm
import DCIPtools as DCIP
from mpl_toolkits.axes_grid1 import make_axes_locatable

# # INPUTS
# LINE/GRID ROOT FOLDER
dataRoot = 'C:\\Users\\HugoLarnier\\Documents\\projects\\orano\\data\\SD0\\'
#dataRoot = 'C:\\Users\\HugoLarnier\\Documents\\projects\\mmg_mcArthur\\L41\\L41\\'
#dataRoot = 'C:\\Users\\HugoLarnier\\Documents\\projects\\india\\3D\\test\\'
#dataRoot = 'C:\\Users\\HugoLarnier\\Documents\\projects\\FieldSchool2018\\FieldSchool\\'

# RECORDINGS FILE
#recordingFile = 'Log\\Recordings_37.txt'
recordingFile = 'Log\\Recordings_41test.txt'
recordingFile = 'Logs\\Recordings_SD0.txt'

# GPS COORDINATES FILES
gpsFile = 'snapCoordinates_CL0.IoBAT.csv'
# gpsFile = 'gpsCoordinates.csv'

# Reading GPS location file
if len(glob(dataRoot + gpsFile)):
    f = open(dataRoot + gpsFile, 'r')
    theoreticalCoordinates = csv.reader(f)
    gpsStations = []
    for row in theoreticalCoordinates:
        try:
            gpsStations.append([float(row[0]), float(row[1])])
        except ValueError:
            pass

# Getting list of nodes in the DATA folder
nodeList = list_nodes(dataRoot + 'DATA\\')

# Checking node library for missing files or empty folders
logFile = open(dataRoot + 'checkDataFolder.log', 'w')
for node in nodeList:
    print("Checking node: " + node)
    logFile.write("Checking node: " + node + '\n')
    datNode = get_list_file(dataRoot + 'DATA', node)
    check_list_files(datNode, logFile)
logFile.close()

# Checking if previous library file is existing
if len(glob(dataRoot + 'nodeLibrary.dat')):
    records = read_library_file(dataRoot + 'nodeLibrary.dat')
    records, nodeList = add_new_nodes(records, nodeList)
else:
    records = [[] for i in range(len(nodeList))]

# Reading recording file
logInjections = read_log_file(dataRoot + recordingFile)

nodeIndex = 0
# Making library of records to compare with injections
for node in nodeList:
    print("Processing node: " + node)
    datNode = get_list_file(dataRoot + 'DATA', node)
    fileIndex = 0
    for datFile in datNode:  # loop on DAT files in folder
        test = is_file_in_library(records, dataRoot + 'DATA\\' + node + '\\' + datFile, nodeIndex)
        if test == 1:
            print("\tReading file: " + datFile)
            fIn = open(dataRoot + 'DATA\\' + node + '\\' + datFile)
            linesFIn = fIn.readlines()
            fIn.close()
            info = get_dat_info(linesFIn)               # Getting DAT file info from header
            timeStamp = get_start_end_time(linesFIn)    # Getting time stamp
            gpsLocations = get_average_gps(linesFIn)    # Getting average GPS location from DAT file
            try:
                records[nodeIndex].append(Record(node_id=info[0],
                                                 node_type=info[3],
                                                 name=dataRoot + 'DATA\\' + node + '\\' + datFile,
                                                 mem=info[1],
                                                 start_date=timeStamp[0],
                                                 end_date=timeStamp[1],
                                                 relay_state=info[2],
                                                 northing=gpsLocations[0],
                                                 easting=gpsLocations[1],
                                                 altitude=gpsLocations[2]))
            except IndexError:
                pass
        else:
            print("\t File: " + datFile + " is already in the library")

        fileIndex += 1
    nodeIndex += 1

# Making updated library file
make_library_file(nodeList, records, dataRoot)

# Comparing recording file and list of records
for injection in logInjections:
    for nodeIndex in range(len(nodeList)):
        for record in records[nodeIndex]:
            test = is_node_active(injection, record)
            if test:
                injection.list_nodes.append(record.node_id)
                injection.list_gps.append([record.northing, record.easting, record.altitude])
                if record.node_type == 'C':
                    injection.list_type.append('C')
                else:
                    injection.list_type.append(record.relay_state)
                injection.list_files.append(record.name)

# Creating synthetic signal for harmonic analysis
print('Starting harmonic analysis of DATA folder')
sample_freq = 150.
sample_half_t = 600.0

# Creating synthetic current to get frequencies of harmonics
time_s, data_s = signal.synthetic_current(sample_freq, 2)

# Getting frequencies of harmonics
periods = signal.get_maxima(1,
                            len(data_s),
                            np.linspace(0, 1, len(data_s)) * sample_freq,
                            np.abs(np.fft.fft(data_s)))

# Small outlier removal
index_2 = signal.get_frequency_index(np.linspace(0, 1, len(data_s)) * sample_freq, periods)
ps_values_2 = signal.get_harmonics_amplitude(index_2,
                                            np.linspace(0, 1, len(data_s)) * sample_freq,
                                            np.abs(np.fft.fft(data_s)))
_, periods = signal.remove_outliers(ps_values_2, periods)

# Parameters for Vp calculation
start_vp = 50
end_vp = 90
sample_freq = 150.
sample_half_t = 600.0

# Doing analysis for each injection
for injection in logInjections:
    start_date = str(injection.start_date.year) + '-' + str(injection.start_date.month) + '-' +\
                 str(injection.start_date.day) + '/' + str(injection.start_date.hour) + ':' +\
                 str(injection.start_date.minute) + ':' + str(injection.start_date.second)
    end_date = str(injection.end_date.year) + '-' + str(injection.end_date.month) + '-' +\
                 str(injection.end_date.day) + '/' + str(injection.end_date.hour) + ':' +\
                 str(injection.end_date.minute) + ':' + str(injection.end_date.second)
    print("Analysis injection " + str(injection.num) + "/" + str(logInjections[-1].num))

    # Initializing arrays for harmonics power and Vp values storage
    values = [np.nan for i in range(len(injection.list_gps))]
    Vp = [np.nan for i in range(len(injection.list_gps))]

    # Switching to UTM for distances calculations
    utm_coord = [[], []]
    for i in range(len(injection.list_gps)):
        tmp = utm.from_latlon(injection.list_gps[i][0], injection.list_gps[i][1])
        utm_coord[0].append(tmp[0])
        utm_coord[1].append(tmp[1])

    # Starting loop on node files for each injection
    for i in range(len(injection.list_gps)):
        test = 1
        print("\tAnalyzing file: " + injection.list_files[i])
        if injection.list_type[i] == 'A':       # Only getting nodes power
            fIn = open(injection.list_files[i], 'r')
            linesFIn = fIn.readlines()
            fIn.close()
            try:
                time, data = read_data(linesFIn)
            except:
                print('\t -> Cannot read file.')
                test = 0

            num_half_T = np.round(np.asarray(data).size / sample_half_t) - 2
            trim_length = num_half_T * sample_half_t
            data = np.asarray(data[0:int(trim_length)])
            time = np.asarray(time[0:int(trim_length)])
            if test:
                try:
                    freq, amp, _ = signal.get_fft(time, data)                    # FFT calculation
                    index = signal.get_frequency_index(freq, periods)            # Get frequency index for harmonics
                    ps_values = signal.get_harmonics_amplitude(index, freq, amp) # Get amplitude of harmonics for DAT file
                except:
                    print('\t Issue with file.')
                    test = 0
            if test:
                values[i] = np.log10(ps_values[0])                               # Switching to log10
            else:
                values[i] = np.nan
            distances = signal.get_distances(utm_coord)                          # Getting distances for future use
            ratios = signal.get_power_ratio(values)                              # Getting power ratio between nodes
            try:
                window = DCIP.createHanningWindow(num_half_T)   # creates filter window
                tHK = DCIP.filterKernal(filtershape=window)     # creates filter kernal
                tHK.kernel = tHK.kernel[1:-1]
                stack = tHK * data                               # stack data
                Vp[i] = DCIP.getPrimaryVoltage(start_vp, end_vp, stack)         # Use of DCIP Tools for VP calculation
            except:
                Vp[i] = np.nan

    name_couples = np.asarray([[injection.list_nodes[i] + "/" + injection.list_nodes[j] for j in range(len(injection.list_nodes))] for i in range(len(injection.list_nodes))]) # For plot purposes
    distances.resize((int(len(utm_coord[0]) ** 2), 1))  # For plot and testing purposes, resizing distances for graph(distance, ratio)
    ratios.resize((int(ratios.shape[0] ** 2), 1))  # For plot and testing purposes, resizing distances for graph(distance, ratio)
    name_couples.resize((int(name_couples.shape[0] ** 2), 1))  # For plot purposes

    # Getting Vp and power max and min values for plotting purposes
    cmin = int(np.floor(np.min([values[i] for i in range(len(values)) if not np.isnan(values[i])])))
    cmax = int(np.ceil(np.max([values[i] for i in range(len(values)) if not np.isnan(values[i])])))
    cminVp = int(np.floor(np.min([Vp[i] for i in range(len(Vp)) if not np.isnan(Vp[i])])))
    cmaxVp = int(np.ceil(np.max([Vp[i] for i in range(len(Vp)) if not np.isnan(Vp[i])])))

    # Plotting
    fig, (ax1, ax2) = plt.subplots(2, 1)
    divider1 = make_axes_locatable(ax1)
    divider2 = make_axes_locatable(ax2)
    fig.suptitle('Injection ' + str(injection.num) + '; Type: ' + injection.valid + '; Date: ' + start_date +
                 ' to ' + end_date)
    # Scatter of harmonic value depending on location to observe spatial variation
    sc1 = ax1.scatter([utm_coord[0][i] for i in range(len(utm_coord[0])) if injection.list_type[i] == 'A'],
                [utm_coord[1][i] for i in range(len(utm_coord[1])) if injection.list_type[i] == 'A'],
                c=[values[i] for i in range(len(injection.list_gps)) if injection.list_type[i] == 'A'], s=100, vmin=cmin, vmax=cmax, cmap=cm.plasma)
    # Scatter of Vp value depending on location to observe spatial variation
    sc2 = ax2.scatter([utm_coord[0][i] for i in range(len(utm_coord[0])) if injection.list_type[i] == 'A'],
                [utm_coord[1][i] for i in range(len(utm_coord[1])) if injection.list_type[i] == 'A'],
                c=[Vp[i] for i in range(len(injection.list_gps)) if injection.list_type[i] == 'A'], s=100, vmin=cminVp, vmax=cmaxVp, cmap=cm.plasma)
    cax1 = divider1.append_axes('right', size='5%', pad=0.05)
    cax2 = divider2.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(sc1, cax=cax1, orientation='vertical')
    fig.colorbar(sc2, cax=cax2, orientation='vertical')
    """
    ax2.plot(distances, ratios, 'ro')
    for i in range(len(distances)):
        ax2.annotate(name_couples[i], (distances[i], ratios[i]))
    """
    for i in range(len(injection.list_nodes)):
        if injection.list_type[i] == 'A':
            #plt.plot(injection.list_gps[i][1], injection.list_gps[i][0], 'k+')
            ax2.plot(utm_coord[0][i], utm_coord[1][i], 'k+')
            ax2.annotate(injection.list_nodes[i], (utm_coord[0][i], utm_coord[1][i]))
            ax1.plot(utm_coord[0][i], utm_coord[1][i], 'k+')
            ax1.annotate(injection.list_nodes[i], (utm_coord[0][i], utm_coord[1][i]))
    plt.show()
exit()


# Plotting injection map
shortLegend = mlines.Line2D([], [], color='magenta', marker='v', linestyle='None', label='Short')
nodeLegend = mlines.Line2D([], [], color='black', marker='o', linestyle='None', label='Node')
currentRecorderLegend = mlines.Line2D([], [], color='red', marker='o', linestyle='None', label='Current Recorder')

figInj, ax0 = plt.subplots(1, 1)
figInj.set_size_inches(8.27, 11.69)
figInj.suptitle('Injection map; Mem numbers from: ' + str(logInjections[0].num) + ' to ' + str(logInjections[-1].num))
for injection in logInjections:
    for i in range(len(injection.list_nodes)):
        if injection.list_type[i] == 'C':
            ax0.plot(injection.list_gps[i][1], injection.list_gps[i][0], 'or')
            ax0.annotate(injection.list_nodes[i], (injection.list_gps[i][1], injection.list_gps[i][0]))

if len(glob(dataRoot + gpsFile)):
    for loc in gpsStations:
        ax0.plot(loc[0], loc[1], 'k+')
    ax0.set_xlabel('Longitude (degrees)')
    ax0.set_ylabel('Latitude (degrees)')
    ax0.grid()
ax0.legend(handles=[currentRecorderLegend], loc='upper center',
           bbox_to_anchor=(0.5, -0.1), fancybox=True, shadow=True, ncol=1)
figInj.tight_layout()
figInj.subplots_adjust(top=0.95)
plt.show()

# Plotting node map for each injection
for injection in logInjections:
    start_date = str(injection.start_date.year) + '-' + str(injection.start_date.month) + '-' +\
                 str(injection.start_date.day) + '/' + str(injection.start_date.hour) + ':' +\
                 str(injection.start_date.minute) + ':' + str(injection.start_date.second)
    end_date = str(injection.end_date.year) + '-' + str(injection.end_date.month) + '-' +\
                 str(injection.end_date.day) + '/' + str(injection.end_date.hour) + ':' +\
                 str(injection.end_date.minute) + ':' + str(injection.end_date.second)
    fig, (ax0, ax1) = plt.subplots(2, 1)
    fig.set_size_inches(8.27, 11.69)
    fig.suptitle('Injection ' + str(injection.num) + '; Type: ' + injection.valid + '; Date: ' + start_date +
                 ' to ' + end_date)
    for i in range(len(injection.list_gps)):
        if injection.list_type[i] == 'A':
            ax0.plot(injection.list_gps[i][1], injection.list_gps[i][0], 'ok')
        elif injection.list_type[i] == 'S':
            ax0.plot(injection.list_gps[i][1], injection.list_gps[i][0], 'v', color='magenta')
        elif injection.list_type[i] == 'C':
            ax0.plot(injection.list_gps[i][1], injection.list_gps[i][0], 'or')
        ax0.annotate(injection.list_nodes[i], (injection.list_gps[i][1], injection.list_gps[i][0]))
    if len(glob(dataRoot + gpsFile)):
        for loc in gpsStations:
            ax0.plot(loc[0], loc[1], 'r+')
    ax0.set_xlabel('Longitude (degrees)')
    ax0.set_ylabel('Latitude (degrees)')
    ax0.grid()
    ax0.legend(handles=[shortLegend, nodeLegend, currentRecorderLegend], loc='upper center',
               bbox_to_anchor=(0.5, -0.1), fancybox=True, shadow=True, ncol=3)
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Current (mA)')
    strTitle = 'Current recorders: '
    count_recorder = 0
    for i in range(len(injection.list_type)):
        if injection.list_type[i] == 'C':
            count_recorder += 1
            strTitle += injection.list_nodes[i] + '; '
            fIn = open(injection.list_files[i], 'r')
            linesFIn = fIn.readlines()
            fIn.close()
            time, data = read_data(linesFIn)
            tmp = injection.list_files[i]
            tmp = tmp.split(dataRoot)
            ax1.plot(time, data - np.mean(data), label=tmp[1])
    ax1.grid()
    ax1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1),
               fancybox=True, shadow=True, ncol=1)
    #ax1.suptitle(strTitle)
    #fig.savefig(dataRoot + 'injection_' + str(injection.num) + '.pdf', papertype='a4', quality=95,
    #            orientation='portrait')
    fig.tight_layout()
    fig.subplots_adjust(top=0.95)
    plt.show()
