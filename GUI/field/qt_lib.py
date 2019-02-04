import sys
from PyQt5 import QtGui, QtWidgets
from PyQt5.QtCore import pyqtSlot
from io_lib import Record, Injection
import io_lib as iolib
import clustering_lib as clustering
from glob import glob
from functools import partial
from matplotlib.pylab import plt
from matplotlib.backends.backend_qt5agg import (
        FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
import matplotlib.lines as mlines
import numpy as np
import utm
import csv
import gpxpy
import gpxpy.gpx

shortLegend = mlines.Line2D([], [], color='magenta', marker='v', linestyle='None', label='Short')
nodeLegend = mlines.Line2D([], [], color='black', marker='o', linestyle='None', label='Node')
currentRecorderLegend = mlines.Line2D([], [], color='red', marker='o', linestyle='None', label='Current Recorder')

class Project:

    def __init__(self,
                 path=None,
                 data_path=None,
                 recording_file=None,
                 list_nodes=None,
                 injections=None,
                 records=None):

        self.path = path
        self.data_path = data_path
        self.recording_file = recording_file
        self.list_nodes = list_nodes
        self.injections = injections
        self.records = records

    def check_path(self):
        tmp = glob.glob(self.path)
        if len(tmp) == 0:
            print("Project path doesn't exist.")
            sys.exit()

        tmp = glob.glob(self.recording_file)
        if len(tmp) == 0:
            print("Recording file doesn't exist.")
            sys.exit()

    def get_node_list(self):
        list_nodes = []
        for i in range(len(self.records)):
            if len(self.records[i]):
                list_nodes.append(self.records[i][0].node_id)
        self.list_nodes = list_nodes

    def initialize_node_data(self):
        self.node_data = [[] for i in range(len(self.injections))]

class locations_plot(QtWidgets.QMainWindow):

    def __init__(self, parent=None):
        super(locations_plot, self).__init__(parent)
        self._main = QtWidgets.QWidget()
        self.title = 'Locations viewer - DIAS GEOPHYSICAL'
        self.left = 10
        self.top = 10
        self.width = 1300
        self.height = 700
        self.project = parent.project
        self.index = 0
        self.initUI()
		
    def loadGPX(self):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        file, _ = QtWidgets.QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*)", options=options)
        if file:
            print('GPX file read:', file)
            try:
                gpx_data = iolib.read_gpx(file)
                self.figure.clf()
                self.ax0 = self.figure.add_subplot(111)
                for info in gpx_data:
                	self.ax0.plot(info[0], info[1], 'r*')
                	self.ax0.annotate(info[2], (info[0], info[1]))
                #self.ax0.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1),
                #           fancybox=True, shadow=True, ncol=1)
						   
                self.ax0.set_ylabel('Northing (km)')
                self.ax0.set_xlabel('Easting (km)')
                self.canvas.draw()
            except:
                print('Issue with the file.')
        else:
            print('Cannot read file.')
		
    def saveGPX(self):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        file, _ = QtWidgets.QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*)", options=options)
        if file:
            print('GPX file read:', file)
        gpx_data = iolib.read_gpx(file)        
	
    def clusterFiles(self):
        gpx = clustering.cluster_data_folder(self.project.recording_file, self.project.data_path)
        directory = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select a folder where to save GPX:', 'C:\\', QtWidgets.QFileDialog.ShowDirsOnly)
        if directory:
            print('Directory read:', directory)
        fOut = open(directory + '/GPX_folder.gpx', 'w')
        fOut.write(gpx.to_xml())
        fOut.close()
        print('File saved.')
		
    def clusterNodesMessages(self):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        file, _ = QtWidgets.QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*)", options=options)
        if file:
            print('File read:', file)
        node_messages_file = file
        gpx = clustering.cluster_nodes_messages(node_messages_file)
        directory = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select a folder where to save GPX:', 'C:\\', QtWidgets.QFileDialog.ShowDirsOnly)
        if directory:
            print('Directory read:', directory)
        fOut = open(directory + '/GPX_nodes_messages.gpx', 'w')
        fOut.write(gpx.to_xml())
        fOut.close()
        print('File saved.')
	    
    def initUI(self):
	    self.index = 0
	    self.setWindowTitle(self.title)
	    self.setWindowIcon(QtGui.QIcon('pythonlogo.png'))
	    self.setGeometry(self.left, self.top, self.width, self.height)

	    self.main_widget = QtWidgets.QWidget(self)
	    self.main_widget.setGeometry(100, 75, self.frameGeometry().width() - 200, self.frameGeometry().height() - 50 * 2)
	    self.layout = QtWidgets.QVBoxLayout(self.main_widget)
	    self.figure = plt.Figure(dpi=80)
	    self.canvas = FigureCanvas(self.figure)
	    self.layout.addWidget(self.canvas)
	    self.addToolBar(NavigationToolbar(self.canvas, self))

	    mainMenu = self.menuBar()
	    fileMenu = mainMenu.addMenu('File')
	    loadSingle = QtWidgets.QAction('Load GPX file', self)
	    loadSingle.setShortcut('Ctrl+O')
	    loadSingle.setStatusTip('Load single GPS file to see track')
	    loadSingle.triggered.connect(self.loadGPX)

	    loadNode = QtWidgets.QAction('Load node folder', self)
	    loadNode.setShortcut('Ctrl+O')
	    loadNode.setStatusTip('Load node folder')
	    #loadNode.triggered.connect(self.plot_node_track)

	    exitButton = QtWidgets.QAction(QtGui.QIcon('exit24.png'), 'Exit', self)
	    exitButton.setShortcut('Ctrl+Q')
	    exitButton.setStatusTip('Exit application')
	    exitButton.triggered.connect(self.close)

	    clusteringMenu = mainMenu.addMenu('Clustering')
	    clusterFiles = QtWidgets.QAction('Cluster node files', self)
	    clusterFiles.setShortcut('Ctrl+N')
	    clusterFiles.setStatusTip('Clustering node DAT files')
	    clusterFiles.triggered.connect(self.clusterFiles)
	    clusterMessages = QtWidgets.QAction('Cluster node messages', self)
	    clusterMessages.setShortcut('Ctrl+M')
	    clusterMessages.setStatusTip('Clustering node messages files')
	    clusterMessages.triggered.connect(self.clusterNodesMessages)		
	    fileMenu.addAction(loadSingle)
	    fileMenu.addAction(loadNode)
	    fileMenu.addAction(exitButton)

	    clusteringMenu.addAction(clusterFiles)
	    clusteringMenu.addAction(clusterMessages)

	    self.show()

class gps_plot(QtWidgets.QMainWindow):

    def __init__(self, parent=None):
        super(gps_plot, self).__init__(parent)
        self._main = QtWidgets.QWidget()
        self.title = 'GPS plot - DIAS GEOPHYSICAL'
        self.left = 10
        self.top = 10
        self.width = 1000
        self.height = 800
        self.project = parent.project
        self.index = 0
        self.initUI()


    def plot_single_track(self):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        file, _ = QtWidgets.QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*);;Python Files (*.py)", options=options)
        if file:
            print('File read:', file)
            fIn = open(file, 'r')
            fLines = fIn.readlines()
            fIn.close()
            infos = iolib.get_gps_constellation(fLines)

            easting = infos[1]
            northing = infos[0]
            gpsAverage = [0., 0.]
            gpsAverage[0] = np.median([easting[i] for i in range(len(easting)) if not np.isnan(easting[i])])
            gpsAverage[1] = np.median([northing[i] for i in range(len(northing)) if not np.isnan(northing[i])])
            tmp = utm.from_latlon(gpsAverage[1], gpsAverage[0])
            gpsAverage[0] = tmp[0]
            gpsAverage[1] = tmp[1]
            for i in range(len(easting)):
                try:
                    tmp = utm.from_latlon(northing[i], easting[i])
                    northing[i] = tmp[1]
                    easting[i] = tmp[0]
                    if np.abs(northing[i] - gpsAverage[1]) > 50 or np.abs(easting[i] - gpsAverage[0]) > 50:
                        easting[i] = np.nan
                        northing[i] = np.nan
                except utm.error.OutOfRangeError:
                    easting[i] = np.nan
                    northing[i] = np.nan
            self.figure.clf()
            self.ax0 = self.figure.add_subplot(311)
            self.ax0.plot(easting, northing, 'ro')
            self.ax0.plot(gpsAverage[0], gpsAverage[1], '*b')
            self.ax0.grid()
            self.ax1 = self.figure.add_subplot(312)
            self.ax1.plot(easting, 'ro')
            self.ax2 = self.figure.add_subplot(313)
            self.ax2.plot(northing, 'bo')
            if len(glob(self.project.gps_file)):
                f = open(self.project.gps_file, 'r')
                theoreticalCoordinates = csv.reader(f)
                gpsStations = [[], []]
                numStation = []
                for row in theoreticalCoordinates:
                    try:
                        numStation.append(str(int(row[0])) + ',' + str(int(row[1])))
                    except:
                        pass
                    try:
                        gpsStations[0].append(float(row[2]))
                        gpsStations[1].append(float(row[3]))
                    except ValueError:
                        pass

                f.close()
                for loc in gpsStations:
                    self.ax0.plot(gpsStations[0], gpsStations[1], 'r+')

            self.canvas.draw()


    def plot_node_track(self):

        options = QtWidgets.QFileDialog.Options()
        options = QtWidgets.QFileDialog.DontUseNativeDialog
        directory = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select a folder:', 'C:\\', QtWidgets.QFileDialog.ShowDirsOnly)
        if directory:
            print('Directory read:', directory)

        dat_list = glob(directory + '/*.DAT')
        self.figure.clf()
        self.ax0 = self.figure.add_subplot(111)
        for datFile in dat_list:
            print("\tProcessing file: " + datFile)
            fIn = open(datFile)
            linesFIn = fIn.readlines()
            fIn.close()
            try:
                gpsAverage = iolib.get_average_gps(linesFIn)
                tmp = utm.from_latlon(gpsAverage[0], gpsAverage[1])
                gpsAverage[0] = tmp[0]
                gpsAverage[1] = tmp[1]
                self.ax0.plot(gpsAverage[0], gpsAverage[1], '*b')
            except:
                pass
        self.ax0.grid()

        if len(glob(self.project.gps_file)):
            f = open(self.project.gps_file, 'r')
            theoreticalCoordinates = csv.reader(f)
            gpsStations = [[], []]
            numStation = []
            for row in theoreticalCoordinates:
                try:
                    numStation.append(str(int(row[0])) + ',' + str(int(row[1])))
                except:
                    pass
                try:
                    gpsStations[0].append(float(row[2]))
                    gpsStations[1].append(float(row[3]))
                except ValueError:
                    pass

            f.close()
            for loc in gpsStations:
                self.ax0.plot(gpsStations[0], gpsStations[1], 'r+')

        self.canvas.draw()

    def initUI(self):
        self.index = 0
        self.setWindowTitle(self.title)
        self.setWindowIcon(QtGui.QIcon('pythonlogo.png'))
        self.setGeometry(self.left, self.top, self.width, self.height)

        self.main_widget = QtWidgets.QWidget(self)
        self.main_widget.setGeometry(250, 100, self.frameGeometry().width() - 300, self.frameGeometry().height() - 100 * 2)
        self.layout = QtWidgets.QVBoxLayout(self.main_widget)
        self.figure = plt.Figure(dpi=80)
        self.canvas = FigureCanvas(self.figure)
        self.layout.addWidget(self.canvas)
        self.addToolBar(NavigationToolbar(self.canvas, self))

        mainMenu = self.menuBar()
        fileMenu = mainMenu.addMenu('File')
        loadSingle = QtWidgets.QAction('Load single file', self)
        loadSingle.setShortcut('Ctrl+O')
        loadSingle.setStatusTip('Load single GPS file to see track')
        loadSingle.triggered.connect(self.plot_single_track)

        loadNode = QtWidgets.QAction('Load node folder', self)
        loadNode.setShortcut('Ctrl+O')
        loadNode.setStatusTip('Load node folder')
        loadNode.triggered.connect(self.plot_node_track)

        exitButton = QtWidgets.QAction(QtGui.QIcon('exit24.png'), 'Exit', self)
        exitButton.setShortcut('Ctrl+Q')
        exitButton.setStatusTip('Exit application')
        exitButton.triggered.connect(self.close)

        fileMenu.addAction(loadSingle)
        fileMenu.addAction(loadNode)
        fileMenu.addAction(exitButton)

        self.show()

class data_plot(QtWidgets.QMainWindow):

    def __init__(self, parent=None):
        super(data_plot, self).__init__(parent)
        self._main = QtWidgets.QWidget()
        self.title = 'GPS plot - DIAS GEOPHYSICAL'
        self.left = 10
        self.top = 10
        self.width = 1000
        self.height = 800
        self.project = parent.project
        self.index = 0
        self.initUI()


    def plot_data(self):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        file, _ = QtWidgets.QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*);;Python Files (*.py)", options=options)
        if file:
            print('File read:', file)
            fIn = open(file, 'r')
            fLines = fIn.readlines()
            fIn.close()
            time, data = iolib.read_data(fLines)
            info = iolib.get_dat_info(fLines)
            self.figure.clf()
            self.figure.suptitle('Mem #: ' + str(info[1]))
            self.ax0 = self.figure.add_subplot(111)
            self.ax0.plot(time, data - np.mean(data), label=file)

            self.ax0.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1),
                   fancybox=True, shadow=True, ncol=1)
            self.canvas.draw()

    def initUI(self):
        self.index = 0
        self.setWindowTitle(self.title)
        self.setWindowIcon(QtGui.QIcon('pythonlogo.png'))
        self.setGeometry(self.left, self.top, self.width, self.height)

        self.main_widget = QtWidgets.QWidget(self)
        self.main_widget.setGeometry(250, 100, self.frameGeometry().width() - 300, self.frameGeometry().height() - 100 * 2)
        self.layout = QtWidgets.QVBoxLayout(self.main_widget)
        self.figure = plt.Figure(dpi=80)
        self.canvas = FigureCanvas(self.figure)
        self.layout.addWidget(self.canvas)
        self.addToolBar(NavigationToolbar(self.canvas, self))

        mainMenu = self.menuBar()
        fileMenu = mainMenu.addMenu('File')
        loadSingle = QtWidgets.QAction('Load single file', self)
        loadSingle.setShortcut('Ctrl+O')
        loadSingle.setStatusTip('Load single data file')
        loadSingle.triggered.connect(self.plot_data)

        exitButton = QtWidgets.QAction(QtGui.QIcon('exit24.png'), 'Exit', self)
        exitButton.setShortcut('Ctrl+Q')
        exitButton.setStatusTip('Exit application')
        exitButton.triggered.connect(self.close)

        fileMenu.addAction(loadSingle)
        fileMenu.addAction(exitButton)

        self.show()




class injection_plot(QtWidgets.QMainWindow):

    def __init__(self, parent=None):
        super(injection_plot, self).__init__(parent)
        self._main = QtWidgets.QWidget()
        self.title = 'Injection plot'
        self.left = 10
        self.top = 10
        self.width = 1000
        self.height = 800
        self.project = parent.project
        self.index = 0
        self.initUI()

    def plot_injection(self, text):
        index = int(text)
        injection = self.project.injections[index]
        start_date = str(injection.start_date.year) + '-' + str(injection.start_date.month) + '-' +\
                     str(injection.start_date.day) + '/' + str(injection.start_date.hour) + ':' +\
                     str(injection.start_date.minute) + ':' + str(injection.start_date.second)
        end_date = str(injection.end_date.year) + '-' + str(injection.end_date.month) + '-' +\
                     str(injection.end_date.day) + '/' + str(injection.end_date.hour) + ':' +\
                     str(injection.end_date.minute) + ':' + str(injection.end_date.second)


        self.figure.clf()
        self.ax0 = self.figure.add_subplot(211)
        self.ax1 = self.figure.add_subplot(212)
        self.figure.suptitle('Injection ' + str(injection.num) + '; Type: ' + injection.valid + '; Date: ' + start_date +
                     ' to ' + end_date)
        for i in range(len(injection.list_gps)):
            if injection.list_type[i] == 'A':
                self.ax0.plot(injection.list_gps[i][0], injection.list_gps[i][1], 'ok')
            elif injection.list_type[i] == 'S':
                self.ax0.plot(injection.list_gps[i][0], injection.list_gps[i][1], 'v', color='magenta')
            elif injection.list_type[i] == 'C':
                self.ax0.plot(injection.list_gps[i][0], injection.list_gps[i][1], 'or')
            self.ax0.annotate(injection.list_nodes[i], (injection.list_gps[i][0], injection.list_gps[i][1]))

        if len(glob(self.project.gps_file)):
            f = open(self.project.gps_file, 'r')
            theoreticalCoordinates = csv.reader(f)
            gpsStations = [[], []]
            numStation = []
            for row in theoreticalCoordinates:
                try:
                    numStation.append(str(int(row[0])) + ',' + str(int(row[1])))
                except:
                    pass
                try:
                    gpsStations[0].append(float(row[2]))
                    gpsStations[1].append(float(row[3]))
                except ValueError:
                    pass

            f.close()
            for loc in gpsStations:
                self.ax0.plot(gpsStations[0], gpsStations[1], 'r+')

        self.ax0.set_xlabel('Easting (m)')
        self.ax0.set_ylabel('Northing (m)')
        self.ax0.grid()
        self.ax0.legend(handles=[shortLegend, nodeLegend, currentRecorderLegend], loc='upper center',
                   bbox_to_anchor=(0.5, -0.1), fancybox=True, shadow=True, ncol=3)
        self.ax1.set_xlabel('Time')
        self.ax1.set_ylabel('Current (mA)')
        strTitle = 'Current recorders: '
        count_recorder = 0
        for i in range(len(injection.list_type)):
            if injection.list_type[i] == 'C':
                count_recorder += 1
                strTitle += injection.list_nodes[i] + '; '
                fIn = open(injection.list_files[i], 'r')
                linesFIn = fIn.readlines()
                fIn.close()
                time, data = iolib.read_data(linesFIn)
                tmp = injection.list_files[i]
                tmp = tmp.split(self.project.data_path)
                self.ax1.plot(time, data - np.mean(data), label=tmp[1])
        self.ax1.grid()
        self.ax1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1),
                   fancybox=True, shadow=True, ncol=1)
        #ax1.suptitle(strTitle)
        #fig.savefig(dataRoot + 'injection_' + str(injection.num) + '.pdf', papertype='a4', quality=95,
        #            orientation='portrait')

        self.figure.tight_layout()
        self.figure.subplots_adjust(top=0.95)

        self.canvas.draw()

    def initUI(self):
        print(len(self.project.injections), self.project.injections[0])
        self.index = 0
        self.setWindowTitle(self.title)
        self.setWindowIcon(QtGui.QIcon('pythonlogo.png'))
        self.setGeometry(self.left, self.top, self.width, self.height)

        self.label = QtWidgets.QLabel('Select injection', self)
        self.label.resize(150, 50)
        self.label.move(50, 100)

        comboBox = QtWidgets.QComboBox(self)
        comboBox.addItem('Select injection')
        for i in range(len(self.project.injections)):
            #comboBox.addItem("Index: " + str(i) + ';MEM: ' + str(self.project.injections[i].num))
            comboBox.addItem(str(i))
        comboBox.resize(200, 50)
        comboBox.move(50, 150)
        comboBox.activated.connect(self.plot_injection)


        self.main_widget = QtWidgets.QWidget(self)
        self.main_widget.setGeometry(250, 100, self.frameGeometry().width() - 300, self.frameGeometry().height() - 100 * 2)
        self.layout = QtWidgets.QVBoxLayout(self.main_widget)
        self.figure = plt.Figure(dpi=80)
        self.canvas = FigureCanvas(self.figure)
        self.layout.addWidget(self.canvas)
        self.addToolBar(NavigationToolbar(self.canvas, self))

        mainMenu = self.menuBar()
        fileMenu = mainMenu.addMenu('File')
        exitButton = QtWidgets.QAction(QtGui.QIcon('exit24.png'), 'Exit', self)
        exitButton.setShortcut('Ctrl+Q')
        exitButton.setStatusTip('Exit application')
        exitButton.triggered.connect(self.close)
        fileMenu.addAction(exitButton)


        self.show()


class App(QtWidgets.QMainWindow):

    def __init__(self):
        super().__init__()
        self._main = QtWidgets.QWidget()
        self.title = 'DIAS Field GUI'
        self.left = 10
        self.top = 10
        self.width = 640 * 2
        self.height = 400 * 2
        self.project = Project()
        self.project.data_path = "C:/Users/HugoLarnier/Desktop/Projects/MMG_McArthur_2018/L36/DATA/"
        self.project.recording_file = "C:/Users/HugoLarnier/Desktop/Projects/MMG_McArthur_2018/L36/LogFiles/Recordings_36.txt"
        self.project.gps_file = "C:/Users/HugoLarnier/Desktop/Projects/MMG_McArthur_2018/L36/LogFiles/Recordings_36.txt"
        self.project.list_nodes = []
        self.project.injections = []

        self.initUI()


    def openDirectoryDialog(self):

        options = QtWidgets.QFileDialog.Options()
        options = QtWidgets.QFileDialog.DontUseNativeDialog
        directory = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select a folder:', 'C:\\', QtWidgets.QFileDialog.ShowDirsOnly)
        if directory:
            print('Directory read:', directory)
        self.project.data_path = directory + '/'

    def openFileDialog(self):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        file, _ = QtWidgets.QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*);;Python Files (*.py)", options=options)
        if file:
            print('File read:', file)
        self.project.recording_file = file

    def openGPSFile(self):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        file, _ = QtWidgets.QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*);;CSV Files (*.csv)", options=options)
        if file:
            print('File read:', file)
        self.project.gps_file = file

    def readFiles(self):
        print(self.project.data_path, self.project.recording_file)
        if self.project.data_path and self.project.recording_file:
            print("Reading...")
            logFile = open(self.project.data_path + 'checkDataFolder.log', 'w')

            self.project.list_nodes = iolib.list_nodes(self.project.data_path)

            for node in self.project.list_nodes:
                print("Checking node: " + node)
                logFile.write("Checking node: " + node + '\n')
                datNode = iolib.get_list_files(self.project.data_path, node)
                iolib.check_list_files(datNode, logFile)
            logFile.close()
            print("Done! You can now check the file " + self.project.data_path + "checkDataFolder.log")

            if len(glob(self.project.data_path + 'nodeLibrary.dat')):
                self.project.records = iolib.read_library_file(self.project.data_path + 'nodeLibrary.dat')
                self.project.records, nodeList = iolib.add_new_nodes(self.project.records, self.project.list_nodes)
            else:
                self.project.records = [[] for i in range(len(self.project.list_nodes))]

            # Reading recording file
            self.project.injections = iolib.read_log_file(self.project.recording_file)

            nodeIndex = 0
            # Making library of records to compare with injections
            for node in self.project.list_nodes:
                print("Processing node: " + node)
                datNode = iolib.get_list_files(self.project.data_path, node)
                fileIndex = 0
                for datFile in datNode:  # loop on DAT files in folder
                    test = iolib.is_file_in_library(self.project.records, self.project.data_path + node + '/' + datFile, nodeIndex)
                    if test == 1:
                        print("\tReading file: " + datFile)
                        fIn = open(self.project.data_path + node + '/' + datFile)
                        linesFIn = fIn.readlines()
                        fIn.close()
                        info = iolib.get_dat_info(linesFIn)               # Getting DAT file info from header
                        timeStamp = iolib.get_start_end_time(linesFIn)    # Getting time stamp
                        gpsLocations = iolib.get_average_gps(linesFIn)    # Getting average GPS location from DAT file
                        try:
                            self.project.records[nodeIndex].append(Record(node_id=info[0],
                                                             node_type=info[3],
                                                             name=self.project.data_path + node + '/' + datFile,
                                                             mem=info[1],
                                                             start_date=timeStamp[0],
                                                             end_date=timeStamp[1],
                                                             relay_state=info[2],
                                                             northing=gpsLocations[0],
                                                             easting=gpsLocations[1],
                                                             altitude=gpsLocations[2],
															 line=info[6],
															 station=info[7]))
                        except IndexError:
                            pass
                    else:
                        print("\t File: " + datFile + " is already in the library")

                    fileIndex += 1
                nodeIndex += 1

            # Making updated library file
            iolib.make_library_file(self.project.list_nodes, self.project.records, self.project.data_path)


            # Comparing recording file and list of records
            for injection in self.project.injections:
                for nodeIndex in range(len(self.project.list_nodes)):
                    for record in self.project.records[nodeIndex]:
                        test = iolib.is_node_active(injection, record)
                        if test:
                            injection.list_nodes.append(record.node_id)
                            tmp = record.getUtm()
                            injection.list_gps.append([tmp[0], tmp[1], record.altitude])
                            if record.node_type == 'C':
                                injection.list_type.append('C')
                            else:
                                injection.list_type.append(record.relay_state)
                            injection.list_files.append(record.name)
        else:
            print("You forgot to load stuff. Please do.")

		
    def plotInjections(self):
        #figInj, ax0 = plt.subplots(1, 1)
        self.figure.clf()
        self.main_ax = self.figure.add_subplot(111)
        #figInj.set_size_inches(8.27, 11.69)
        self.figure.suptitle('Injection map; Mem numbers from: ' + str(self.project.injections[0].num) + ' to ' + str(self.project.injections[-1].num))
        for injection in self.project.injections:
            for i in range(len(injection.list_nodes)):
                if injection.list_type[i] == 'C':
                    self.main_ax.plot(injection.list_gps[i][0], injection.list_gps[i][1], 'or')
                    self.main_ax.annotate(injection.list_nodes[i], (injection.list_gps[i][0], injection.list_gps[i][1]))
        if len(glob(self.project.gps_file)):
            f = open(self.project.gps_file, 'r')
            theoreticalCoordinates = csv.reader(f)
            gpsStations = [[], []]
            numStation = []
            for row in theoreticalCoordinates:
                try:
                    numStation.append(str(int(row[0])) + ',' + str(int(row[1])))
                except:
                    pass
                try:
                    gpsStations[0].append(float(row[2]))
                    gpsStations[1].append(float(row[3]))
                except ValueError:
                    pass

            f.close()
            for loc in gpsStations:
                self.main_ax.plot(gpsStations[0], gpsStations[1], 'r+')


        self.main_ax.grid()
        self.main_ax.legend(handles=[currentRecorderLegend], loc='upper center',
                   bbox_to_anchor=(0.5, -0.1), fancybox=True, shadow=True, ncol=1)
        self.figure.tight_layout()
        self.figure.subplots_adjust(top=0.95)


        self.canvas.draw()

    def plotInjection_single(self, static_canvas):

        ex = injection_plot(self)

    def gps_plot_window(self):

        ex2 = gps_plot(self)

    def data_plot_window(self):
    	
    	    ex3 = data_plot(self)

    def clusterData(self):
	    ex4 = locations_plot(self)
		
    def checkMissingData(self):
	    options = QtWidgets.QFileDialog.Options()
	    options |= QtWidgets.QFileDialog.DontUseNativeDialog
	    file, _ = QtWidgets.QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*);;CSV Files (*.csv)", options=options)
	    if file:
	        print('File read:', file)
	    iolib.checkMissingFiles(file, self.project.data_path)


    def initUI(self):
        self.setWindowTitle(self.title)
        self.setWindowIcon(QtGui.QIcon('pythonlogo.png'))
        self.setGeometry(self.left, self.top, self.width, self.height)

        mainMenu = self.menuBar()
        fileMenu = mainMenu.addMenu('File')
        viewMenu = mainMenu.addMenu('View')
        toolsMenu = mainMenu.addMenu('Tools')
        helpMenu = mainMenu.addMenu('Help')

        loadButton = QtWidgets.QAction(QtGui.QIcon('exit24.png'), 'Load DATA', self)
        loadButton.setShortcut('Ctrl+O')
        loadButton.setStatusTip('Load DATA folder')
        loadButton.triggered.connect(self.openDirectoryDialog)

        recordingButton = QtWidgets.QAction(QtGui.QIcon('exit24.png'), 'Load recording file', self)
        recordingButton.setShortcut('Ctrl+O')
        recordingButton.setStatusTip('Load DATA folder')
        recordingButton.triggered.connect(self.openFileDialog)

        loadGPSFile = QtWidgets.QAction('Read GPS file', self)
        loadGPSFile.setShortcut('Ctrl+G')
        loadGPSFile.setStatusTip('Load GPS file for display')
        loadGPSFile.triggered.connect(self.openGPSFile)

        readButton = QtWidgets.QAction('Read folder', self)
        readButton.setShortcut('Ctrl+R')
        readButton.setStatusTip('Read DATA folder')
        readButton.triggered.connect(self.readFiles)

        exitButton = QtWidgets.QAction(QtGui.QIcon('exit24.png'), 'Exit', self)
        exitButton.setShortcut('Ctrl+Q')
        exitButton.setStatusTip('Exit application')
        exitButton.triggered.connect(self.close)

        fileMenu.addAction(loadButton)
        fileMenu.addAction(recordingButton)
        fileMenu.addAction(readButton)
        fileMenu.addAction(loadGPSFile)
        fileMenu.addAction(exitButton)


        self.main_widget = QtWidgets.QWidget(self)
        self.main_widget.setGeometry(50, 100, self.frameGeometry().width() - 50 * 2, self.frameGeometry().height() - 100 * 2)
        self.layout = QtWidgets.QVBoxLayout(self.main_widget)
        self.figure = plt.Figure(dpi=80)
        self.canvas = FigureCanvas(self.figure)
        self.layout.addWidget(self.canvas)
        self.addToolBar(NavigationToolbar(self.canvas, self))

        injectionButton = QtWidgets.QAction('Injection map', self)
        injectionButton.setStatusTip('Plot injection map')
        injectionButton.triggered.connect(partial(self.plotInjections))

        recordButton = QtWidgets.QAction('Single injection map', self)
        recordButton.setStatusTip('Plot injection map')
        recordButton.triggered.connect(self.plotInjection_single)

        gpsButton = QtWidgets.QAction('GPS information', self)
        gpsButton.setStatusTip('Plot GPS information')
        gpsButton.triggered.connect(self.gps_plot_window)

        dataButton = QtWidgets.QAction('Data information', self)
        dataButton.setStatusTip('Data GPS information')
        dataButton.triggered.connect(self.data_plot_window)

        viewMenu.addAction(injectionButton)
        viewMenu.addAction(recordButton)
        viewMenu.addAction(gpsButton)
        viewMenu.addAction(dataButton)

        gpsClusteringButton = QtWidgets.QAction('GPS locations', self)
        gpsClusteringButton.setStatusTip('GPS locations')
        gpsClusteringButton.triggered.connect(partial(self.clusterData))

        missingDataButton = QtWidgets.QAction('Check missing data', self)
        missingDataButton.setStatusTip('Check for missing data')
        missingDataButton.triggered.connect(partial(self.checkMissingData))

        toolsMenu.addAction(gpsClusteringButton)
        toolsMenu.addAction(missingDataButton)
        self.show()
