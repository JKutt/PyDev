import numpy as np
import DCIPtools as DCIP

# ============================ User area ==============================
# ver 0.1
# set variables
start_time = 40
end_time = 1920
direction = 'x'
a = 50
projectName = "Boggy_Creek"
date = "22/03/2019"
line = "8800"
fname = "E:\\Projects\\debug\\Charlotte\\Boggy_L8800_Prelim_QC.DAT"
outname = "E:\\Projects\\debug\\Charlotte\\Boggy_L8800_Prelim_QC-Geosoft.DAT"

# load the data file
patch = DCIP.loadDias(fname)

# write the file
patch.writeOldFartDat(outname, start_time, end_time, direction, a, projectName, date, line)

print("done!!!! damn these 2D old school people!!!")