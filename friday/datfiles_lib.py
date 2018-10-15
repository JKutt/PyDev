## IMPORTS
import numpy as np
import datetime
from os.path import isdir, getsize
from glob import glob
from sys import exit

## Written by H. LARNIER, Oct 2018 - DIAS GEOPHYSICAL LTD.

class Injection:
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
                 altitude=None):
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


def get_list_file(path, node_id):
    # List all dat files for the node "node_id" in the folder "path"
    # Inputs:
    # path: string containing the path to the node to investigate
    # node_id: string containing the node id (two characters)
    # Outputs:
    # dat_list: list of DAT files

    dat_list = glob(path + '\\' + node_id + '\\*.DAT')
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
        tmp = dat_list[i].split('\\')
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
        # Recovering numbers
        for file in list_files:
            tmp = file.split('.')
            num.append(int(tmp[0][2:]))

        diff = [num[i + 1] - num[i] for i in range(len(num) - 1)]
        if np.max(diff) > 1:
            tmp = [i for i in range(len(diff)) if diff[i] > 1]
            for i in range(len(tmp)):
                print('\tWarning: Possible missing files between: ' + str(list_files[tmp[i]]) +
                      ' and ' + str(list_files[tmp[i] + 1]))
                log.write('\tWarning: Possible missing files between: ' + str(list_files[tmp[i]]) +
                          ' and ' + str(list_files[tmp[i] + 1]) + '\n')
        else:
            pass
    else:
        print('\tWarning: Empty folder')
        log.write('\tWarning: Empty folder' + '\n')
    return


def get_dat_info(lines):
    # Get information of DAT file from header
    # Input:
    # lines: list of lines from DAT file
    # Output:
    # dat_info: list of information (Unit ID, MEM number, Relay state, if current of potential)
    dat_info = ['', 0, '', '']
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
    for line in lines:
        if line[0:6] == '$GNGGA' or line[0:6] == '$GPGGA':
            [test_n, test_w, test_a] = [0, 0, 0]
            tmp = line.split(',')
            try:
                north.append(gps_from_minutes_to_decimal(float(tmp[2]) / 100.))
                if tmp[3] == 'S':
                    north[-1] *= -1
            except ValueError:
                test_n = 1
                pass
            except IndexError:
                test_n = 1
                pass
            try:
                west.append(gps_from_minutes_to_decimal(float(tmp[4]) / 100.))
                if tmp[5] == 'W':
                    west[-1] *= -1
            except ValueError:
                test_w = 1
                pass
            except IndexError:
                test_w = 1
                pass
            try:
                alt.append(float(tmp[9]))
            except ValueError:
                test_a = 1
                pass
            except IndexError:
                test_a = 1
                pass
            if test_n and len(north):
                del north[-1]
            if test_w and len(west):
                del west[-1]
            if test_a and len(alt):
                del alt[-1]

    if not len(north):
        return []
    else:
        return [north, west, alt]


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
            tmp = line.split(',')
            tmp = tmp[1]
            if len(tmp) == 9 or len(tmp) == 10:
                try:
                    hour_tmp = int(tmp[0:2])
                except ValueError:
                    print(tmp)
                    exit()
                try:
                    minute_tmp = int(tmp[2:4])
                except ValueError:
                    print(tmp)
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


def list_nodes(data_folder):
    # Get list of nodes from DATA folder
    # Input:
    # data_folder: String containing the DATA folder location
    # Output:
    # node_list: list of nodes (have to be only two characters long)
    node_list = glob(data_folder + '??')
    not_nodes = []
    for line in node_list:
        if not isdir(line):
            not_nodes.append(line)

    for i in range(len(not_nodes)):
        node_list.remove(not_nodes)

    for i in range(len(node_list)):
        node_list[i] = node_list[i][-2:]

    return node_list


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


def read_data(lines):
    # Read data from DAT file
    # Input:
    # lines: list of lines from the DAT file
    # Output:
    # data: list of data
    data = []
    factor = 1
    for line in lines:
        if line[1:11] == 'Conversion':
            tmp = line.split(':')
            factor = int(tmp[1][:-1])
        if line[0] == '+' or line[0] == '-':
            data.append(int(line[:-1]))

    for i in range(len(data)):
        data[i] /= factor
    return data


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
            node_list.append(line[1:3])

    return node_list


def read_library_file(library_file):

    file = open(library_file, 'r')
    lines = file.readlines()
    node_list = get_node_list_from_library(lines)
    recs = [[] for i in range(len(node_list))]
    count = 0
    for line in lines:
        if line[0] == '<':  # new node
            node_id = line[1:3]
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
            count += 1
        tmp = None
        if count == 9:
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
                                               altitude=float(altitude)))
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


def is_file_in_library(records, dat_file, node_index):

    test = 1
    for i in range(len(records[node_index])):
        if dat_file == records[node_index][i].name:
            test = 0
            break

    return test
