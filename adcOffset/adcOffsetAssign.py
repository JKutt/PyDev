from os import listdir
from os.path import isfile, join, isdir
from multiprocessing import Pool

# ==================================================================
#  Global variables
csv_path = "C:/Users/johnk\Projects/orano-nodebiases.csv"
dir_path = "E:/Orano/IPDB/DATA/"
# ==================================================================


def processNode(path):
    #  ===============================================================
    #  read csv and make dictionary

    high_gain_dict = {}
    low_gain_dict = {}
    csv = open(csv_path)
    node_id = []
    high_gain = []
    low_gain = []
    for dline in csv:
        splitt = dline.split(",")
        node_id.append(splitt[0])
        high_gain.append(splitt[5])
        low_gain.append(splitt[6])
    csv.close()
    # print(low_gain)
    for i in range(len(node_id)):
        high_gain_dict[node_id[i]] = high_gain[i]
        low_gain_dict[node_id[i]] = low_gain[i]

    tmp_dir_path = dir_path + path
    onlyfiles = [f for f in listdir(tmp_dir_path)
                 if isfile(join(tmp_dir_path, f))]     # gets files
    for i in range(len(onlyfiles)):
        xyz = open(tmp_dir_path + "/" + onlyfiles[i])  # open to read
        x = []
        for line in xyz:
            x.append((line))
        xyz.close()
        # write IP data to file
        extension = onlyfiles[i].split(".")           # determine extension
        if extension[1:][0] == 'DAT':
            out_file = open(tmp_dir_path +
                            "/" + onlyfiles[i], "w")  # open for writing
            temp_id = "  "
            conversion = 0.
            offset = 0.
            for lines in range(len(x)):
                peices = x[lines].split()             # look for adc
                if len(peices) > 1:
                    if peices[0] == "#Unit":
                        temp_id = peices[2]
                    if peices[0] == "#Conversion":
                        conversion_ = peices[2].split(":")
                        conversion = float(conversion_[1])
                    if peices[0] == "#ADC":
                        try:
                            if conversion > 300000.:
                                offset = high_gain_dict[temp_id]
                            else:
                                offset = low_gain_dict[temp_id]
                        except KeyError:
                            pass
                            # print("could not match ID")
                        x[lines] = "#ADC offset: " + str(offset) + "\n"
                out_file.write("%s" % x[lines])        # write to file
            out_file.close()                          # close file


if __name__ == '__main__':
    #  ===================================================================
    #  begin adc offset assign
    print("assigning offsets")
    only_nodes = [f for f in listdir(dir_path) if isdir(join(dir_path, f))]
    pool = Pool()
    pool.map(processNode, only_nodes)
    print("done.... go home!")
