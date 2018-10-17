#!C:\\Users\\HugoLarnier\\Anaconda3\\python

from datfiles_lib import *
import csv
import matplotlib.pylab as plt
import matplotlib.lines as mlines
from glob import glob

# # INPUTS
# LINE/GRID ROOT FOLDER
dataRoot = 'C:\\Users\\HugoLarnier\\Documents\\projects\\orano\\data\\CL0E\\'
#dataRoot = 'C:\\Users\\HugoLarnier\\Documents\\projects\\india\\3D\\test\\'
#dataRoot = 'C:\\Users\\HugoLarnier\\Documents\\projects\\FieldSchool2018\\FieldSchool\\'

# RECORDINGS FILE
recordingFile = 'LOGS\\Recordings_0.txt'

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

# Plotting injection map
shortLegend = mlines.Line2D([], [], color='magenta', marker='v', linestyle='None', label='Short')
nodeLegend = mlines.Line2D([], [], color='black', marker='o', linestyle='None', label='Node')
currentRecorderLegend = mlines.Line2D([], [], color='red', marker='o', linestyle='None', label='Current Recorder')
"""
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
"""
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
            ax1.plot(time, data, label=tmp[1])
    ax1.grid()
    ax1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1),
               fancybox=True, shadow=True, ncol=1)
    #ax1.suptitle(strTitle)
    #fig.savefig(dataRoot + 'injection_' + str(injection.num) + '.pdf', papertype='a4', quality=95,
    #            orientation='portrait')
    fig.tight_layout()
    fig.subplots_adjust(top=0.95)
    plt.show()