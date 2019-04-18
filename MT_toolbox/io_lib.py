import os
import numpy as np
import sys
import array


def write_to_log_file(log_file, message):
    f_out = open(log_file, 'a')
    f_out.write(message)
    f_out.close()

def read_binary_file(file_in):
    f_in = open(file_in, 'rb')
    s = f_in.read()
    file_out = array.array('d', s)
    file_out = np.asarray(file_out)

    return file_out

def read_ascii_file(file_in):
    file_out = []
    f_in = open(file_in, 'r')
    lines = f_in.readlines()
    for line in lines:
        file_out.append(float(line[:-1]))

    f_in.close()
    return file_out

def read_phoenix_file(file_in):

    Ex = []
    Ey = []
    Hx = []
    Hy = []
    Hz = []
    f_in = open(file_in, 'r')
    lines = f_in.readlines()
    for line in lines:
        if "Number of channels:" in line:
            print(line.split("Number of channels:")[1][:-1])
            n_channels = int(line.split("Number of channels:")[1][:-1])
        if "Scans per record:" in line:
            sample_freq = int(line.split("Scans per record:")[1][:-1])

        try:
            n = int(line[:11])
            Ex.append(int(line[:21]))
            Ey.append(int(line[:32]))
            Hx.append(int(line[:43]))
            Hy.append(int(line[:54]))
            Hz.append(int(line[:65]))
        except:
            pass
    f_in.close()
    print(Ex)
    return Ex, Ey, Hx, Hy, Hz
