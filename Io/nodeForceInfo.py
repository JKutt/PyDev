from os import listdir
from os.path import isfile, join, isdir
# from multiprocessing import Pool

# ==================================================================
#  Setting Variables
dir_path = ["/Users/juan/testData/G2000263.DAT"]  # use this for specific

# USE THIS FOR BULK NODE DIRECTORY PROCESSING
# node_dir = "C:/Users/johnk/Projects/DATA/VR/"
# dir_path = [f for f in listdir(node_dir)
#                    if isfile(join(node_dir, f))]     # get files

shift = None                          # set amount of shift
lat = None                             # set forced latitude
longi = None                           # set forced longitude
elev = None                            # sets elevation
# ==================================================================


def processNode(path,
                option=None,
                timeshift=None,
                longitude=None,
                latitude=None,
                elevation=None):
    #  ===============================================================
    #  read file to check if force headers already exisit
    flag_shift = False                               # flag that shift present
    flag_position = False                            # flag that pos. present
    flag_elevation = False                           # flags elevation
    for index in range(len(path)):
        csv = open(path[index])
        for dline in csv:
            splitt = dline.split(":")
            if splitt[0] == "#Time Shift":
                flag_shift = True
            if splitt[0] == "#Override Easting":
                flag_position = True
            if splitt[0] == "#Override Elevation":
                flag_elevation = True

        csv.close()
        if flag_shift:
            print("overwrite force time")           # mention overwrite
        if flag_position:
            print("overwrite force pos.")           # mention overwrite

        xyz = open(path[index])  # open to read
        x = []
        for line in xyz:
            x.append((line))
        xyz.close()
        print("size of x: ", len(x))
        # write desired force data to file
        cnt = 0
        extension = path[index].split(".")           # determine extension
        if extension[1:][0] == 'DAT':
            out_file = open(path[index], "w")            # open for writing
            if timeshift is not None:
                if flag_shift:                        # shift present
                    for lines in range(len(x)):
                        peices = x[lines].split(":")      # look for adc
                        if len(peices) > 1:
                            if peices[0] == "#Time Shift":
                                x[lines] = "#Time Shift:" + str(timeshift) + "\n"
                        out_file.write("%s" % x[lines])   # write to file

                else:
                    force_time = "#Time Shift:" + str(timeshift) + "\n"
                    out_file.write("%s" % force_time)
                    for lines in range(len(x)):
                        out_file.write("%s" % x[lines])
            if longitude is not None:
                if flag_position:                        # shift present
                    for lines in range(len(x)):
                        peices = x[lines].split(":")      # look for adc
                        if len(peices) > 1:
                            if peices[0] == "#Override Easting":
                                x[lines] = "#Override Easting:" + str(longitude) + "\n"
                            if peices[0] == "#Override Northing":
                                x[lines] = "#Override Northing:" + str(latitude) + "\n"
                        out_file.write("%s" % x[lines])   # write to file

                else:
                    force_long = "#Override Easting:" + str(latitude) + "\n"
                    force_lat = "#Override Northing:" + str(longitude) + "\n"
                    out_file.write("%s" % force_long)
                    out_file.write("%s" % force_lat)
                    for lines in range(len(x)):
                        out_file.write("%s" % x[lines])

            if elevation is not None:
                if flag_elevation:                        # shift present
                    for lines in range(len(x)):
                        peices = x[lines].split(":")      # look for adc
                        if len(peices) > 1:
                            if peices[0] == "#Override Elevation":
                                x[lines] = "#Override Elevation:" + str(elevation) + "\n"
                        out_file.write("%s" % x[lines])   # write to file

                else:
                    force_elev = "#Override Elevation:" + str(elevation) + "\n"
                    out_file.write("%s" % force_elev)
                    for lines in range(len(x)):
                        out_file.write("%s" % x[lines])

            out_file.close()                          # close file


if __name__ == '__main__':
    #  ===================================================================
    #  begin adc offset assign
    print("looking at raw file")
    processNode(dir_path,
                timeshift=shift,
                latitude=lat,
                longitude=longi,
                elevation=elev)
