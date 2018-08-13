from os import listdir
from os.path import isfile, join, isdir
from multiprocessing import Pool
import numpy as np

# written by JK Aug 2018 V0.1
# ==================================================================
#  Global variables
csv_path = "S:/Projects/Aurora/Malva/IPDB/11425/Locations_L11425_iobat.csv"
dir_path = "S:/Projects/Aurora/Malva/IPDB/11425/DATA2/"
# ==================================================================

def processNode(path):
    #  ===============================================================
    #  read csv and make dictionary
    node_dict = {}
    csv = open(csv_path)
    node_id = []
    latitude = []
    longitude = []
    easting = []
    northing = []
    for dline in csv:
        splitt = dline.split(",")
        splitt2 = splitt[7].split("\n")
        node_id.append(splitt2[0])
        latitude.append(splitt[5])
        longitude.append(splitt[6])
    csv.close()
    # print(node_id)

    for i in range(len(node_id)):
        node_dict[node_id[i]] = i

    # print(node_dict["XR"])
    tmp_dir_path = dir_path + path
    onlyfiles = [f for f in listdir(tmp_dir_path)
                 if isfile(join(tmp_dir_path, f))]     # gets files
    for i in range(len(onlyfiles)):
        extension = onlyfiles[i].split(".")           # determine extension
        if extension[1:][0] == 'DAT':
            xyz = open(tmp_dir_path + "/" + onlyfiles[i])  # open to read
            x = []
            flag_position = False                            # flag that pos. present
            lat = []
            lon = []
            for line in xyz:
                x.append((line))
                splitt = line.split(":")
                splitt2 = line.split(",")
                if splitt[0] == "#Override Easting":
                    flag_position = True
                if splitt2[0] == "#GNRMC":
                    try:
                        lat = float(splitt2[3]);
                        flat = floor(lat / 100.);
                        minutes = ((lat / 100.) - flat) * 100.;
                        declat = flat + minutes / 60.;
                        # check the hemishpere
                        if splitt2[4] == "S":
                            declat = -1 * declat;
                        lat.append(declat)
                    except:
                        pass
                    try:
                        lon = float(splitt2[5]);
                        flat = floor(lon / 100.);
                        minutes = ((lon / 100.) - flat) * 100.;
                        declon = flat + minutes / 60.;
                        # check the hemishpere
                        if splitt2[6] == "W":
                            declon = -1 * declon;
                        lon.append(declon)
                    except:
                        pass

            xyz.close()
            # write IP data to file
            lat = np.asarray(lat)
            lon = np.asarray(lon)
            lat_file = np.median(lat)
            lon_file = np.median(lon)
            
            id_found = False
            node_index = 0
            curr_dist = 1000000.
            for jdex in range(len(latitude)):
                dist = np.sqrt((latitude[jdex] - lat_file)**2 +
                 (longitude[jdex] - lon_file)**2) 
                if dist < curr_dist:
                    file_id = onlyfiles[i][:2]
                    if file_id == node_id[jdex]
                        node_index = jdex
                        curr_dist = dist
                        id_found = True

            if id_found: 
                out_file = open(tmp_dir_path +
                                "/" + onlyfiles[i], "w")  # open for writing
                if flag_position:                        # shift present
                    for lines in range(len(x)):
                        peices = x[lines].split(":")      # look for adc
                        if len(peices) > 1:
                            if peices[0] == "#Override Easting":
                                x[lines] = "#Override Easting:" + longitude[node_index] + "\n"
                            if peices[0] == "#Override Northing":
                                x[lines] = "#Override Northing:" + latitude[node_index] + "\n"
                        out_file.write("%s" % x[lines])   # write to file

                else:
                    force_long = "#Override Easting:" + latitude[node_index] + "\n"
                    force_lat = "#Override Northing:" + longitude[node_index] + "\n"
                    out_file.write("%s" % force_long)
                    out_file.write("%s" % force_lat)
                    for lines in range(len(x)):
                        out_file.write("%s" % x[lines])
                out_file.close()                          # close file
            


if __name__ == '__main__':
    #  ===================================================================
    #  begin adc offset assign
    print("assigning offsets")
    only_nodes = [f for f in listdir(dir_path) if isdir(join(dir_path, f))]
    pool = Pool()
    pool.map(processNode, only_nodes)
    # for i in range(len(only_nodes)):
    #     processNode(only_nodes[i])
    print("done.... go home!")
