from os import listdir
from os.path import isfile, join, isdir
from multiprocessing import Pool

# ==================================================================
#  Global variables
csv_path = "S:/Projects/Aurora/Malva/IPDB/11425/Locations_L11425_iobat.csv"
dir_path = "S:/Projects/Aurora/Malva/IPDB/11425/DATA/"
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
            for line in xyz:
                x.append((line))
                splitt = line.split(":")
                if splitt[0] == "#Override Easting":
                    flag_position = True
            xyz.close()
            # write IP data to file
            
            id_found = True
            node_index = 0
            node_index = (node_dict[onlyfiles[i][:2]])
            try:
                node_index = (node_dict[onlyfiles[i][:2]])
            except KeyError:
                # print("hit shorted")
                id_found = False
                pass
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
