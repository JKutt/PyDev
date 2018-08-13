import struct
import os
import numpy as np

def getDipoleStack(path):
        bin_file = path

        f = open(bin_file, "rb")          # open the file
        lines = f.read()                  # read data as string (can modify this to get header info
        return_char = 0.                  # beginning of binary data idex
        for idx in range(100):
            if lines[idx] == "\r\n":
                return_char = idx
        bin_start = 178 
        # return_char + 2       # increment 2 spaces to start of binary data
        print(len(lines))
        data = []                         # initiate array for data
        with open(bin_file, 'rb') as f:   # open file for data extraction
            f.seek(int(bin_start),
                   os.SEEK_SET)          # seek to beginning of binary data
            while True:
                b = f.read(8)             # read 8 bytes at a time
                if not b:                 # break out if end of file
                    # eof
                    break
                data.append(
                    struct.unpack('d',
                                  b))  # store data
                if len(b) < 8:
                    break
        return data


file_path = "C:/Projects/Aurora/Volo/IPDB/Grid3/DATA/TXR/DF004585.raw"  
in_data = getDipoleStack(file_path)
in_data = np.asarray(in_data)

print("In: {0}".format(np.sum(in_data[150:290])))