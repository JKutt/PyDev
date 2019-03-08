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
