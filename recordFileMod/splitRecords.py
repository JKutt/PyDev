import numpy as np
from datetime import datetime
import time


class recordInfo:

    def __init__(self, recordnumber, start):
        self.record_number = recordnumber
        self.start_time = start

    def setStop(self, stop):
        self.end_time = stop

    def setStatus(self, status):
        self.status = status

    def getSerialStartTime(self):
        split_start = self.start_time.split("T")
        #  split the date
        split_date = split_start[0].split("-")
        # split time
        split_time = split_start[1].split("Z")
        # lets get some seconds out of this
        string_date = split_date[0] + "/" + split_date[1] + "/" + split_date[2]
        string_date_time = string_date + " " + split_time[0]
        d = datetime.strptime(string_date_time, "%Y/%m/%d %H:%M:%S")
        time_in_sec = time.mktime(d.timetuple())
        return time_in_sec

    def getSerialStopTime(self):
        split_start = self.end_time.split("T")
        #  split the date
        split_date = split_start[0].split("-")
        # split time
        split_time = split_start[1].split("Z")
        # lets get some seconds out of this
        string_date = split_date[0] + "/" + split_date[1] + "/" + split_date[2]
        string_date_time = string_date + " " + split_time[0]
        d = datetime.strptime(string_date_time, "%Y/%m/%d %H:%M:%S")
        time_in_sec = time.mktime(d.timetuple())
        return time_in_sec

    def makeSubReadings(self, numberofreadings):
        start_of_rec = self.getSerialStartTime()
        end_of_rec = self.getSerialStopTime()
        number_of_sub_readings = np.floor((end_of_rec -
                                          start_of_rec) / numberofreadings)
        new_starts = []
        for i in range(int(numberofreadings)):
            rr = start_of_rec + (i * number_of_sub_readings)
            new_starts.append(datetime.fromtimestamp(rr).strftime("%Y-%m-%dT%H:%M:%SZ"))
        return new_starts

# ================================================================
# doing stuff here
rec_file_path = "C:/Users/johnk/devProjects/Python/recordFileMod/Recordings_1800.txt"
rec_file = open(rec_file_path)
recs = []
cnt = -1
current_record = -1
for line in rec_file:
    splitt = line.split()
    if current_record != float(splitt[1]):
        # print("now here first")
        cnt = cnt + 1
        current_record = float(splitt[1])
        record_read = recordInfo(float(splitt[1]), splitt[3])
        recs.append(record_read)
    else:
        if splitt[0] == "Stop:":
            recs[cnt].setStop(splitt[3])
        if splitt[0] == "Status:":
            recs[cnt].setStatus(splitt[3])

# ff = recs[0].getSerialStartTime()
# rr = recs[0].getSerialStopTime()
# ll = datetime.fromtimestamp(rr).strftime("%Y-%m-%dT%H:%M:%SZ")
print(recs[0].makeSubReadings(6.0))
