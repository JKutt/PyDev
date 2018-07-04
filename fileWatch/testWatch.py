import sys
from PyQt5 import QtCore


def directory_changed(path):
    print('Directory Changed: %s' % path)


def file_changed(path):
    print('File Changed: %s' % path)


app = QtCore.QCoreApplication(sys.argv)

paths = ['/Users/juan/Documents/PyDev/fileWatch/']

fs_watcher = QtCore.QFileSystemWatcher(paths)
fs_watcher.directoryChanged.connect(directory_changed)
fs_watcher.fileChanged.connect(file_changed)

sys.exit(app.exec_())
