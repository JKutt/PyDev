# Lib for input/output for friday

from classes import Project, Node
from os import mkdir
from glob import glob
import shutil
import multiprocessing
from sys import exit

def create_quarantine_directory(project):

    if not len(glob(project.path + '\\data_quarantine')):
        try:
            mkdir(project.path + '\\data_quarantine')
        except OSError as exc:
            if exc.errno == errno.ENOSPC:
                raise
                print("No space left.")
                exit()
            elif exc.errno == errno.EEXIST:
                raise
                print("Directory already exists.")


def create_node_directory(project, node):

    if not len(glob(project.path + '\\data_quarantine\\' + node)):
        try:
            mkdir(project.path + '\\data_quarantine\\' + node)
        except OSError as exc:
            if exc.errno == errno.ENOSPC:
                raise
                print("No space left.")
                exit()
            elif exc.errno == errno.EEXIST:
                raise
                print("Directory already exists.")


def quarantine_file(project, node):

    messages = []
    # test that directory exists
    name = node.file_name[-12:]
    if len(glob(project.path + '\\data_quarantine\\' + node.id)) == 0:
        print('Directory does not exist, creating')
        create_node_directory(project, node.id)
    else:
        print('Directory exist')

    print('[' + multiprocessing.current_process().name + '] \t Moving file: ' + name + ' to quarantine.')
    messages.append('[' + multiprocessing.current_process().name + '] \t Moving file: ' + name + 'to quarantine.\n')
    print('-> Moving from: ' + node.file_name + ' to ' + project.path + '\\data_quarantine\\' + node.id + '\\' + name)
    messages.append('[' + multiprocessing.current_process().name + '] \t -> Moving from: ' + node.file_name + ' to ' + project.path + '\\data_quarantine\\' + node.id + '\\' + name + '.\n')
    shutil.move(node.file_name, project.path + '\\data_quarantine\\' + node.id + '\\' + name)

    return messages
