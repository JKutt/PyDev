## IMPORTS
import numpy as np
import datetime
from os.path import isdir, getsize
from glob import glob
from sys import exit
import matplotlib.pylab as plt
import math
import utm
from sys import platform
import julian

import gpxpy
import gpxpy.gpx
from gpxpy.parser import GPXParser

## Written by H. LARNIER, Oct 2018 - DIAS GEOPHYSICAL LTD.

if platform == 'linux1' or platform == 'linux2' or platform == 'darwin':
    separator = '/'
elif platform == 'win32':
    separator = '\\'
else:
    print('Unknown platform')


class Record:
    def __init__(self,
                 name=None,
                 node_id=None,
                 node_type=None,
                 mem=None,
                 start_date=None,
                 end_date=None,
                 relay_state=None,
                 northing=None,
                 easting=None,
                 altitude=None,
				 line=None,
				 station=None):
        self.name = name
        self.node_id = node_id
        self.node_type = node_type
        self.mem = mem
        self.start_date = start_date
        self.end_date = end_date
        self.relay_state = relay_state
        self.northing = northing
        self.easting = easting
        self.altitude = altitude
        self.line = line
        self.station = station

    def getUtm(self):
        tmp = utm.from_latlon(self.northing, self.easting)
        return tmp[0], tmp[1]




class Injection:
    """ List attributes of an injection.

    attributes:
    valid: Noise, Good or Bad
    num: Mem number
    start_date/end_date
    list_nodes: All nodes active during that injection
    list_files: All files related to that injection
    list_gps: GPS coordinates for all files related to that injection
    list_type: Node type for all files related to that injection

    """

    def __init__(self,
                 valid=None,
                 num=None,
                 start_date=None,
                 end_date=None,
                 list_nodes=None,
                 list_files=None,
                 list_gps=None,
                 list_type=None):
        self.valid = valid
        self.num = num
        self.start_date = start_date
        self.end_date = end_date
        self.list_nodes = list_nodes
        self.list_files = list_files
        self.list_gps = list_gps
        self.list_type = list_type



def list_nodes(data_folder):
    # Get list of nodes from DATA folder
    # Input:
    # data_folder: String containing the DATA folder location
    # Output:
    # node_list: list of nodes (have to be only two characters long)
    node_list = glob(data_folder + '*')

    not_nodes = []
    for line in node_list:
        if not isdir(line):
            not_nodes.append(line)
    for i in range(len(not_nodes)):
        node_list.remove(not_nodes[i])

    for i in range(len(node_list)):
        tmp = node_list[i].split(separator)
        node_list[i] = tmp[-1]

    return node_list


def read_gpx(file):	
	xml = open(file, 'r')
	gpx = GPXParser()
	gpx.init(xml)
	gpx.parse()
	infos = []
	for waypoint in gpx.gpx.waypoints:
		try:
			tmp = utm.from_latlon(waypoint.latitude, waypoint.longitude)
			infos.append([tmp[0], tmp[1], waypoint.name])
		except:
			print('Point: ' + waypoint.name + 'cannot be read')
	return infos

def read_log_file(log_file):
    # Reads the log file and return injection number + time stamp
    # Input:
    # log_file: Name of the log file
    # Output:
    # injection: list of injection (class defined above)
    f_in = open(log_file, 'r')
    lines = f_in.readlines()
    f_in.close()
    injections = []
    for i in range(len(lines)):
        if lines[i][0:5] == 'Start':
            tmp = lines[i].split(':')
            num1 = tmp[1]
            date = tmp[2]
            date = date.split('T')
            tmp1 = date[0].split('-')
            tmp2 = date[1].split(':')
            date_start = datetime.datetime(int(tmp1[0]), int(tmp1[1]), int(tmp1[2]),
                                           int(tmp2[0]), int(tmp[3]), int(tmp[4][0:2]))
            tmp = lines[i + 1].split(':')
            num2 = tmp[1]
            date = tmp[2]
            date = date.split('T')
            tmp1 = date[0].split('-')
            tmp2 = date[1].split(':')
            date_end = datetime.datetime(int(tmp1[0]), int(tmp1[1]), int(tmp1[2]),
                                         int(tmp2[0]), int(tmp[3]), int(tmp[4][0:2]))
            tmp = lines[i + 2].split(':')
            num = tmp[1]
            if num == num1 and num == num2:
                valid = tmp[2][:-1]
                injections.append(Injection(valid=valid,
                                            num=int(num),
                                            start_date=date_start,
                                            end_date=date_end,
                                            list_files=[],
                                            list_nodes=[],
                                            list_gps=[],
                                            list_type=[]))
            else:
                injections.append(Injection(valid='Bad',
                                            num=int(num),
                                            start_date=date_start,
                                            end_date=date_end,
                                            list_files=[],
                                            list_nodes=[],
                                            list_gps=[],
                                            list_type=[]))

    return injections

def get_list_files(path, node_id):
    # List all dat files for the node "node_id" in the folder "path"
    # Inputs:
    # path: string containing the path to the node to investigate
    # node_id: string containing the node id (two characters)
    # Outputs:
    # dat_list: list of DAT files


    dat_list = glob(path + '/' + node_id + '/' + '*.DAT')
    to_remove = []
    # Getting empty files
    for i in range(len(dat_list)):
        size = getsize(dat_list[i])
        if size == 0:
            to_remove.append(dat_list[i])

    # Removing empty files
    for i in range(len(to_remove)):
        dat_list.remove(to_remove[i])

    # Only keeping the file name without the path
    for i in range(len(dat_list)):
        tmp = dat_list[i].split(separator)
        dat_list[i] = tmp[-1]

    return dat_list


def check_list_files(list_files, log):
    # Simple checks on the list of files
    # Input:
    # list_files: List of files in node folder
    # log: pointer on log file to write stuff inside
    # Outputs:
    # None
    num = []
    if len(list_files) > 1:
        ## Comparing node_messages_list and file list
        #if len(node_messages_list) > 0:
        #    print()

        # Recovering numbers
        for file in list_files:
            tmp = file.split('.')
            num.append(int(tmp[0][2:]))

        diff = [num[i + 1] - num[i] for i in range(len(num) - 1)]
        if np.max(diff) > 1:
            tmp = [i for i in range(len(diff)) if diff[i] > 1]
            for i in range(len(tmp)):
                log.write('\tWarning: Possible missing files between: ' + str(list_files[tmp[i]]) +
                          ' and ' + str(list_files[tmp[i] + 1]) + '\n')
        else:
            pass
    else:
        log.write('\tWarning: Empty folder' + '\n')
    return


def get_dat_info(lines):
    # Get information of DAT file from header
    # Input:
    # lines: list of lines from DAT file
    # Output:
    # dat_info: list of information (Unit ID, MEM number, Relay state, if current of potential)
    dat_info = ['', 0, '', '', 0., 0., 0, 0]
    for line in lines:
        if not line[0] == '#':
            break
        else:
            if line[1:5] == 'Unit':
                tmp = line.split(':')
                dat_info[0] = tmp[1][:-1].replace(" ", "")
            if line[1:4] == 'Mem':
                tmp = line.split(':')
                dat_info[1] = int(tmp[1][:-1])
            if line[1:6] == 'Relay':
                tmp = line.split(':')
                relay_state = tmp[1][:-1].replace(" ", "")
                dat_info[2] = relay_state
            if line[1:8] == 'Current':
                dat_info[3] = 'C'
            if line[1:8] == 'Voltage':
                dat_info[3] = 'V'
            if line[1:17] == 'Override Easting':
                tmp = line.split(':')
                dat_info[5] = float(tmp[1][:-1])
            if line[1:18] == 'Override Northing':
                tmp = line.split(':')
                dat_info[4] = float(tmp[1][:-1])
            if line[1:5] == 'Line':
                tmp = line.split(':')
                dat_info[6] = int(tmp[1][:-1])
            if line[1:8] == 'Station':
                tmp = line.split(':')
                dat_info[7] = int(tmp[1][:-1])
    return dat_info


def get_gps_constellation(lines):
    # Read the GPS constellation from a DAT file
    # Inputs:
    # dat_file: list of lines from read .DAT file
    # Outputs:
    # list of north, west and elevation

    alt = []
    north = []
    west = []
    quality = []
    number_satellites = []
    infos = [[] for i in range(5)]
    for line in lines:
        if line[0:6] == '$GNGGA' or line[0:6] == '$GPGGA':
            [test_n, test_w, test_a, test_quality, test_number_of_satellites] = [0, 0, 0, 0, 0]
            tmp = line.split(',')
            try:
                infos[0].append(gps_from_minutes_to_decimal(float(tmp[2]) / 100.))
                if tmp[3] == 'S':
                    infos[0][-1] *= -1
            except ValueError:
                test_n = 1
                pass
            except IndexError:
                test_n = 1
                pass
            try:
                infos[1].append(gps_from_minutes_to_decimal(float(tmp[4]) / 100.))
                if tmp[5] == 'W':
                    infos[1][-1] *= -1
            except ValueError:
                test_w = 1
                pass
            except IndexError:
                test_w = 1
                pass
            try:
                infos[2].append(float(tmp[9]))
            except ValueError:
                test_a = 1
                pass
            except IndexError:
                test_a = 1
                pass
            # GPS quality
            try:
                infos[3].append(int(tmp[6]))
            except ValueError:
                test_quality = 1
            except IndexError:
                test_quality = 1

            try:
                infos[4].append(int(tmp[7]))
            except ValueError:
                test_number_of_satellites = 1
            except IndexError:
                test_number_of_satellites = 1

            indexes = [test_n, test_w, test_a, test_quality, test_number_of_satellites]

            if 1 in indexes:
                for i in range(len(infos)):
                    if indexes[i] == 0 and len(infos[i]):
                        del infos[i][-1]
            #print(len(infos[0]), len(infos[1]), len(infos[2]), len(infos[3]) ,len(infos[4]), indexes)


    if not len(infos[0]):
        return []
    else:
        return infos

def get_average_gps(lines):
    # Return the median gps location from the constellation
    # Inputs:
    # dat_file: .DAT file to read
    # Outputs:
    # list of average gps [north, west, elev]
    gps_constellation = get_gps_constellation(lines)
    if len(gps_constellation):
        north = np.median(gps_constellation[0])
        west = np.median(gps_constellation[1])
        elev = np.median(gps_constellation[2])
        return [north, west, elev]
    else:
        return []

def get_time(lines):
    # Read the time from a DAT file
    # Inputs:
    # dat_file: .DAT file read
    # Outputs:
    # list of north, west and elevation

    [year, month, day, hour, minute, sec] = [[], [], [], [], [], []]
    time = []
    for line in lines:
        if line[0:5] == '#Date':
            tmp = line.index(':')
            tmp = line[tmp+1:].split('-')
            year = int(tmp[0])
            month = int(tmp[1])
            day = int(tmp[2])
        if line[0:5] == '#Time':
            tmp = line.index(':')
            tmp = line[tmp + 1:].split(':')
            hour = int(tmp[0])
            minute = int(tmp[1])
            sec = int(tmp[2])
        if not len(time):
            try:
                time.append(datetime.datetime(year, month, day, hour, minute, sec))
            except TypeError:
                pass
            except ValueError:
                time = []
                break


        if line[0:7] == '$GNRMC,' or line[0:7] == '$GPRMC,':
            time_tmp = get_date_from_gps_value(line)
            if not time_tmp == datetime.datetime(2000, 1, 1, 1, 1, 1):
                time.append(time_tmp)
            """
            tmp = line.split(',')
            tmp2 = tmp[]
            tmp = tmp[1]
            if len(tmp) == 9 or len(tmp) == 10:
                try:
                    hour_tmp = int(tmp[0:2])
                except ValueError:
                    exit()
                try:
                    minute_tmp = int(tmp[2:4])
                except ValueError:
                    exit()
                try:
                    sec_tmp = int(tmp[4:6])
                except ValueError:
                    exit()
                try:
                    time.append(datetime.datetime(year, month, day, hour_tmp, minute_tmp, sec_tmp))
                except TypeError:
                    pass
                except ValueError:
                    pass
                """
    return time


def get_start_end_time(lines):
    # Return the end and start date from a DAT file
    # Input:
    # lines: lines from a DAT file
    # Output:
    # List of end and start time from the DAT file
    time = get_time(lines)
    try:
        start_end = [time[0], time[-1]]
        return start_end
    except:
        print('Error with the .DAT file.')
        return []

def write_node_entry(rec, file):
    # Description goes here

    start_date = str(rec.start_date.year) + '-' + str(rec.start_date.month) + '-' + str(rec.start_date.day) \
        + 'T' + str(rec.start_date.hour) + ':' + str(rec.start_date.minute) + ':' + str(rec.start_date.second)
    end_date = str(rec.end_date.year) + '-' + str(rec.end_date.month) + '-' + str(rec.end_date.day) \
        + 'T' + str(rec.end_date.hour) + ':' + str(rec.end_date.minute) + ':' + str(rec.end_date.second)
    file.write('\t<file_name:' + rec.name + '>\n')
    file.write('\t\t<node_type:' + rec.node_type + '>\n')
    file.write('\t\t<node_mem:' + str(rec.mem) + '>\n')
    file.write('\t\t<start_date:' + start_date + '>\n')
    file.write('\t\t<end_date:' + end_date + '>\n')
    file.write('\t\t<relay_state:' + rec.relay_state + '>\n')
    file.write('\t\t<northing:' + str(rec.northing) + '>\n')
    file.write('\t\t<easting:' + str(rec.easting) + '>\n')
    file.write('\t\t<altitude:' + str(rec.altitude) + '>\n')
    file.write('\t\t<line:' + str(rec.line) + '>\n')
    file.write('\t\t<station:' + str(rec.station) + '>\n')

    return


def make_library_file(node_list, recs, root):
    # Make library file so nothing has to be redone again

    node_library = open(root + 'nodeLibrary.dat', 'w')
    for i in range(len(node_list)):
        if len(recs[i]):
            node_library.write('<' + node_list[i] + '>\n')
            for rec in recs[i]:
                write_node_entry(rec, node_library)
    node_library.close()

    return


def get_date_from_library(date):
    # Description goes here
    date = date.split('T')
    tmp1 = date[0].split('-')
    tmp2 = date[1].split(':')
    date_datetime = datetime.datetime(int(tmp1[0]), int(tmp1[1]), int(tmp1[2]),
                                      int(tmp2[0]), int(tmp2[1]), int(tmp2[2]))

    return date_datetime


def get_node_list_from_library(lines):

    node_list = []
    for line in lines:
        if line[0] == '<':  # new node
            node_list.append(line[1:-2])

    return node_list


def read_library_file(library_file):

    file = open(library_file, 'r')
    lines = file.readlines()
    node_list = get_node_list_from_library(lines)
    recs = [[] for i in range(len(node_list))]
    count = 0
    for line in lines:
        if line[0] == '<':  # new node
            node_id = line[1:-2]
            tmp = [i for i in range(len(node_list)) if node_id == node_list[i]]
            node_index = tmp[0]
        else:  # same node, new file, \t,
            if count == 0:
                tmp = line.split('file_name:')
                file_name = tmp[1][:-2]
            elif count == 1:
                tmp = line.split('node_type:')
                node_type = tmp[1][:-2]
            elif count == 2:
                tmp = line.split('node_mem:')
                node_mem = tmp[1][:-2]
            elif count == 3:
                tmp = line.split('start_date:')
                start_date = get_date_from_library(tmp[1][:-2])
            elif count == 4:
                tmp = line.split('end_date:')
                end_date = get_date_from_library(tmp[1][:-2])
            elif count == 5:
                tmp = line.split('relay_state:')
                relay_state = tmp[1][:-2]
            elif count == 6:
                tmp = line.split('northing:')
                northing = tmp[1][:-2]
            elif count == 7:
                tmp = line.split('easting:')
                easting = tmp[1][:-2]
            elif count == 8:
                tmp = line.split('altitude:')
                altitude = tmp[1][:-2]
            elif count == 9:
                tmp = line.split('line:')
                line_data = tmp[1][:-2]
            elif count == 10:
                tmp = line.split('station:')
                station = tmp[1][:-2]
            count += 1
        if count == 11:
            try:
                recs[node_index].append(Record(node_id=node_id,
                                               node_type=node_type,
                                               name=file_name,
                                               mem=int(node_mem),
                                               start_date=start_date,
                                               end_date=end_date,
                                               relay_state=relay_state,
                                               northing=float(northing),
                                               easting=float(easting),
                                               altitude=float(altitude),
											   line=int(line_data),
											   station=int(station)))
            except IndexError:
                pass
            count = 0
    file.close()

    return recs


def get_node_list_from_records(records):

    node_list = []
    for i in range(len(records)):
        if len(records[i]):
            node_list.append(records[i][0].node_id)

    return node_list


def add_new_nodes(records, node_list):

    node_list_library = get_node_list_from_records(records)

    for i in range(len(node_list)):
        j = 0
        while node_list[i] > node_list_library[j]:
            j += 1
        if not node_list[i] == node_list_library[j]:
            node_list_library.insert(j, node_list[i])
            records.insert(j, [])

    return records, node_list_library



def gps_from_minutes_to_decimal(gps_value):
    # Get decimal coordinates from decimal minutes
    # Inputs:
    # gps_value: gps coordinates in decimal minutes
    # Output:
    # gps value in decimal value
    degree = np.floor(gps_value)
    tmp = (gps_value - degree) * 100
    minutes = np.floor(tmp)
    seconds = ((tmp - minutes) * 60.)

    return degree + minutes / 60. + seconds / 3600.


def is_file_in_library(records, dat_file, node_index):

    test = 1
    for i in range(len(records[node_index])):
        if dat_file == records[node_index][i].name:
            test = 0
            break

    return test



def get_date_from_gps_value(line):
    tmp = line.split(',')
    time = datetime.datetime(2000, 1, 1, 1, 1, 1)
    if len(tmp) >= 10:
        tmp2 = tmp[9]
        tmp = tmp[1]
        try:
            year = int("20" + tmp2[4:])
            month = int(tmp2[2:4])
            day = int(tmp2[0:2])
        except ValueError:
            return time
        try:
            hour_tmp = int(tmp[0:2])
        except ValueError:
            return time
        try:
            minute_tmp = int(tmp[2:4])
        except ValueError:
            return time
        try:
            sec_tmp = int(tmp[4:6])
        except ValueError:
            return time
        try:
            time = datetime.datetime(year, month, day, hour_tmp, minute_tmp, sec_tmp)
            return time
        except TypeError:
            return time
        except ValueError:
            return time
    else:
        return time


def is_node_active(inj, rec):
    # Return boolean indicating if the considered node is active during the considered injection
    # Inputs:
    # inj: injection item (class Injection)
    # rec: recording item (class Record)
    # Output:
    # is_active: boolean value indicating if node active during recording or not
    if inj.start_date > rec.end_date or inj.end_date < rec.start_date:
        is_active = False
    else:
        is_active = True

    return is_active

def get_date_from_gps_value(line):
    tmp = line.split(',')
    time = datetime.datetime(2000, 1, 1, 1, 1, 1)
    if len(tmp) >= 10:
        tmp2 = tmp[9]
        tmp = tmp[1]
        try:
            year = int("20" + tmp2[4:])
            month = int(tmp2[2:4])
            day = int(tmp2[0:2])
        except ValueError:
            return time
        try:
            hour_tmp = int(tmp[0:2])
        except ValueError:
            return time
        try:
            minute_tmp = int(tmp[2:4])
        except ValueError:
            return time
        try:
            sec_tmp = int(tmp[4:6])
        except ValueError:
            return time
        try:
            time = datetime.datetime(year, month, day, hour_tmp, minute_tmp, sec_tmp)
            return time
        except TypeError:
            return time
        except ValueError:
            return time
    else:
        return time


def read_data(lines):
    # Read data from DAT file
    # Input:
    # lines: list of lines from the DAT file
    # Output:
    # data: list of data

    data = []           # list of data values
    time_pps = []           # list of time values assigned to data list
    factor = 1

    count = 0
    adc_offset = 0.
    have_timing = 0
    gps_count = 0
    threshold = 20000000.
    assigned_pps = 0
    time_shift = 0.0
    start_time = -1
    last_gps_time = 0.
    have_s = 0.
    l_pps = 0.
    sample_rate = 20.
    gps_count = 0.
    f_pps = 0
    hold_data = 0
    val = 0
    lsval = 0
    for line in lines:
        if line[1:4] == "Time":
            tmp = line.split(':')
            time_shift = (int(tmp[1][:-1]) / 1000.0) / (60. * 60. *24.)
        if line[1:11] == 'Conversion':
            tmp = line.split(':')
            factor = int(tmp[1][:-1])
        if line[1:4] == 'ADC':
            tmp = line.split(':')
            try:
                adc_offset = float(tmp[1][:-1])
            except:
                adc_offset = 0.
        if line[1:7] == 'Sample':
            tmp = line.split(':')
            sample_rate = int(tmp[1][:-1])              # Frequency sample rate
            t_inc = (1. / sample_rate) / (60. * 60. * 24.)  # Time increment between two samples (in years)
        if line[:-1] == 'PPS':
            assigned_pps = 1
            #print("Found PPS")
        if line[0] == 'S':
            #print("Found second")
            tmp = line.split('S')
            val = float(tmp[1][:-1])
            if val < 360000000:
                hold_data = 0
                have_timing += 1
                if have_timing == 0 or gps_count == 0:
                    f_pps = count
                    l_pps = count
                    lsval = count
                else:
                    obs_samplerate = count - l_pps
                    have_s = 1
                    #print("Sample rates", obs_samplerate, sample_rate)
                    if not obs_samplerate == sample_rate:
                        spacing = val - lsval
                        if spacing > 0:
                            samples = count - l_pps                 # Number of samples
                            exp_samples = spacing * sample_rate     # Number of samples from file
                            missing = int(exp_samples - samples)         # Difference is missing samples
                            if missing > 0:
                                tmp = data[-1]
                                for i in range(missing):
                                    data.append(tmp)
                                    count += 1
                            else:
                                count += missing

                    l_pps = count
                    lsval = val
            else:
                hold_data = 1
        if line[0] == '$':
            tmp = line.split(',')
            if tmp[0] == "$GNRMC" or tmp[0] == "$GPRMC" or tmp[0] == "$PSRFC" \
                    or tmp[0] == "$PSRMC" or tmp[0] == "$PPRMC" :
                #print('Reading GPS String')
                gps_valid = 0
                if have_timing > 2 and assigned_pps == 1:
                    assigned_pps = 0
                    if len(tmp) >= 10:  # Check if line is complete
                        try:
                            time_tmp = get_date_from_gps_value(line)
                            time_tmp = julian.to_jd(time_tmp)#get_julian_datetime(time_tmp)
                            #time_tmp = get_julian_datetime(time_tmp)
                            time_tmp += ((time_shift / 1000.) / 60. * 1 / 60. * 1 / 24.)
                            gps_count += 1
                        except ValueError:
                            pass
                        gps_valid = 1
                    if gps_valid == 1:
                        if start_time == -1:
                            start_time = time_tmp - (t_inc * f_pps)
                            time_pps.append(start_time)
                        else:
                            time_pos = l_pps + 1
                            if time_pos < len(data) and have_s == 1:
                                if time_tmp == last_gps_time:
                                    time_tmp = time_tmp + 1. / 86400.
                                if len(time_pps) < l_pps:
                                    while len(time_pps) < l_pps:
                                        time_pps.append(0.)
                                    time_pps[l_pps - 1] = time_tmp
                                else:
                                    time_pps[l_pps - 1] = time_tmp
                                have_s = 0
                                last_gps_time = time_tmp
        if line[0] == '+' or line[0] == '-' or line[0] == ' ':
            if have_timing > 0:
                if hold_data == 0:
                    if len(line) <= 12:
                        try:
                            if int(line[:-1]) / factor - adc_offset < threshold:
                                count += 1
                                data.append(int(line[:-1]))
                        except:
                            pass

    s_time = start_time
    tmp = start_time
    for i in range(len(time_pps)):
        if not time_pps[i] == 0:
            if tmp > time_pps[i]:
                time_pps[i] = 0.

    if len(data) > 0:
        if len(time_pps):
            for i in range(len(time_pps)):
                if not time_pps[i] == 0:
                    s_time = time_pps[i]
                else:
                    s_time = s_time + t_inc
                    time_pps[i] = s_time

                #s_time = s_time + t_inc
                #time_pps[i] = s_time

    if not time_shift == 0:
        for i in range(len(time_pps)):
            time_pps[i] += time_shift

    for i in range(len(data)):
        data[i] = (data[i] / factor) - adc_offset

    return time_pps[0:], data[0:len(time_pps)]

def checkMissingFiles(nodes_messages_file, data_folder):

	fIn = open(nodes_messages_file, 'r')
	lines = fIn.readlines()
	fIn.close()

	list_files = []

	for line in lines:
		if '.dat' in line:
			tmp = line.split(',')
			for string in tmp:
				if '.dat' in string:
					list_files.append(string[:-4])
					
			
	listNodes = list_nodes(data_folder)
	list_upload = []

	for node in listNodes:
		files = get_list_files(data_folder, node)
		for file in files:
			list_upload.append(file[:-4])
	missing = []
	for file in list_files:
		if not file in list_upload:
			missing.append(file)
	missing = np.unique(missing)
	missing = missing.tolist()

	missingNodesList = []
	missingNodesFiles = []
	for file in missing:
		id = file[0:2]
		if not id in missingNodesList:
			missingNodesList.append(id)
			missingNodesFiles.append([])
	missingNodesList = np.sort(missingNodesList)
	missingNodesList = missingNodesList.tolist()
	for file in missing:
		id = file[0:2]
		ind = missingNodesList.index(id)
		missingNodesFiles[ind].append(file + '.dat')
		
		
	fOut = open('list_missing.csv', 'w')
	for i in range(len(missingNodesList)):
		fOut.write(missingNodesList[i] + ',')
	fOut.write('\n')
	nb = [len(missingNodesFiles[i]) for i in range(len(missingNodesList))]
	nb2 = np.max(nb)
	for i in range(nb2):
		for j in range(len(missingNodesList)):
			if i == 1 and nb[j] == 1:
				fOut.write(',')
			else:
				if i >= nb[j]:
					fOut.write(',')
				else:
					fOut.write(missingNodesFiles[j][i] + ',')
		fOut.write('\n')
	print('You can now check the list of missing files at:' + 'list_missing.csv')