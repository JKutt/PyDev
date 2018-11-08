import sys
from PyQt5 import QtGui, QtWidgets
from PyQt5.QtCore import pyqtSlot
from io_lib import Record, Injection
import io_lib as iolib
from glob import glob
from functools import partial
from matplotlib.pylab import plt
from matplotlib.backends.backend_qt5agg import (
        FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
import matplotlib.lines as mlines
import numpy as np

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


class injection_plot(QtWidgets.QMainWindow):

    def __init__(self, parent=None):
        super(injection_plot, self).__init__(parent)
        self._main = QtWidgets.QWidget()
        self.title = 'Injection plot'
        self.left = 10
        self.top = 10
        self.width = 2000
        self.height = 1600
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
        """
        if len(glob(self.project.data_path + gpsFile)):
            for loc in gpsStations:
                ax0.plot(loc[0], loc[1], 'r+')
        """
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
        self.figure = plt.Figure()
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
        self.width = 640 * 3
        self.height = 400 * 3
        self.project = Project()
        self.project.data_path = "C:/Users/HugoLarnier/Desktop/Projects/MMG_McArthur_2018/L36/DATA/"
        self.project.recording_file = "C:/Users/HugoLarnier/Desktop/Projects/MMG_McArthur_2018/L36/LogFiles/Recordings_36.txt"
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

    def readFiles(self):
        print(self.project.data_path, self.project.recording_file)
        if self.project.data_path and self.project.recording_file:
            print("Reading...")
            logFile = open(self.project.data_path + 'checkDataFolder.log', 'w')

            self.project.list_nodes = iolib.list_nodes(self.project.data_path)

            for node in self.project.list_nodes:
                print("Checking node: " + node)
                logFile.write("Checking node: " + node + '\n')
                datNode = iolib.get_list_files(self.project.data_path + 'DATA', node)
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
                                                             altitude=gpsLocations[2]))
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

    def printFiles(self):
        print(vars(self.project))

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


        self.main_ax.grid()
        self.main_ax.legend(handles=[currentRecorderLegend], loc='upper center',
                   bbox_to_anchor=(0.5, -0.1), fancybox=True, shadow=True, ncol=1)
        self.figure.tight_layout()
        self.figure.subplots_adjust(top=0.95)


        self.canvas.draw()
        print('test')

    def plotInjection_single(self, static_canvas):

        ex = injection_plot(self)

    def initUI(self):
        self.setWindowTitle(self.title)
        self.setWindowIcon(QtGui.QIcon('pythonlogo.png'))
        self.setGeometry(self.left, self.top, self.width, self.height)

        mainMenu = self.menuBar()
        fileMenu = mainMenu.addMenu('File')
        viewMenu = mainMenu.addMenu('View')
        searchMenu = mainMenu.addMenu('Search')
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
        fileMenu.addAction(exitButton)


        self.main_widget = QtWidgets.QWidget(self)
        self.main_widget.setGeometry(50, 100, self.frameGeometry().width() - 50 * 2, self.frameGeometry().height() - 100 * 2)
        self.layout = QtWidgets.QVBoxLayout(self.main_widget)
        self.figure = plt.Figure()
        self.canvas = FigureCanvas(self.figure)
        self.layout.addWidget(self.canvas)
        self.addToolBar(NavigationToolbar(self.canvas, self))

        injectionButton = QtWidgets.QAction('Injection map', self)
        injectionButton.setStatusTip('Plot injection map')
        injectionButton.triggered.connect(partial(self.plotInjections))

        recordButton = QtWidgets.QAction('Single injection map', self)
        recordButton.setStatusTip('Plot injection map')
        recordButton.triggered.connect(self.plotInjection_single)

        viewMenu.addAction(injectionButton)
        viewMenu.addAction(recordButton)

        self.show()
