from os import listdir
from os.path import isfile, join

dir_path = "C:/Users/johnk/Projects/CobreGrandeL1/L1/DATA/VU/"
onlyfiles = [f for f in listdir(dir_path) if isfile(join(dir_path, f))]
print("starting timeshift")
for i in range(len(onlyfiles)):
    xyz = open(dir_path + onlyfiles[i])
    x = []
    for line in xyz:
        x.append((line))
    xyz.close()
    # write IP data to file
    out_file = open(dir_path + onlyfiles[i], "w")
    for j in range(len(x)):
        out_file.write("%s" % x[j])
        if j == 2:
            out_file.write("%s\n" % "#Timeshift:-26.666667")
    out_file.close()
    print("timeshift done")
