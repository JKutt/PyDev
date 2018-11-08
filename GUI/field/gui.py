import sys
from PyQt5 import QtGui, QtWidgets
from PyQt5.QtCore import pyqtSlot
from qt_lib import App, Project

if __name__ == '__main__':

    app = QtWidgets.QApplication(sys.argv)
    ex = App()

    sys.exit(app.exec_())
