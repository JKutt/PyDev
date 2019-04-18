from sys import exit, path

import io_lib
import spectrum_lib
import matplotlib.pylab as plt
from classes_ts import project, station_mt, step_parameters
import spectrum_lib
from argparse import ArgumentParser
import datetime
import time


version = '0.1'
date_now = datetime.datetime.now()
date_now = str(date_now.year) + '-' + str(date_now.month) + '-' + str(date_now.day)

parser = ArgumentParser(usage='diasMT: Toolbox for MT data processing',
                        description='Version ' + version,
                        epilog='If no argument is input, this help will be displayed.')
parser.add_argument("-p", "--plot",
                    action='store_true', default=False,
                    help="Read files, plot time series, and exist")

parser.add_argument("-c", "--compute",
                    action='store_true', default=False,
                    help="Read files and compute mt tensor")

parser.add_argument("-t", "--test",
                    action='store_true', default=False,
                    help="Read parameters file and tests it.")

args = parser.parse_args()
if not args.test and not args.plot and not args.compute:
    parser.parse_args(['-h'])

if __name__ == '__main__':
    print('DIAS MT Toolbox - Version ' + str(version))
    ## Initializing project class containing processing information
    project = project(version=version,
                      date_creation=date_now)
    project.read_parameters_file()
    project.check_project()

    ## Reading data
    # Initalizing station class
    station = station_mt(output_level=project.output_level,
                                 log_file=project.log_file)
    station.read_parameters_file_and_load_data()
    station.check_station(project.type)
    if project.parameters.wavelets == 'yes':
        station.get_wavelet_transforms(project.parameters)

    if args.test:
        print('[INFO] Test option was selected. Exiting now...')
        exit()

    if args.plot:
        print('[INFO] Now plotting, and then exiting.')
        station.plot()
        exit()

    regression_time = time.time()
    if args.compute:
        if project.parameters.wavelets == 'yes':
            print('[INFO] Wavelets selected. Now detecting EM waves before starting response function computation.')
            station.detect_elf(project.parameters)
            print('[INFO] ' + str(len(station.events[0])) + ' events detected.')
            #print(len(station.events[0]))
            for time_tmp in station.events[0]:
                plt.plot(station.input[0].data[time_tmp - 100:time_tmp+100])
                plt.plot(station.input[1].data[time_tmp - 100:time_tmp+100], 'r')
                plt.show()

        print('[INFO] Now starting computation of response functions.')
        print('       According to the given parameters, I have ' + str(project.parameters.nb_reductions + 1) + ' iterations to perform.')
        station.load_frequencies(project.parameters)
        step_param = step_parameters()
        for step in range(project.parameters.nb_reductions + 1):
            step_param.update_param(project.parameters, step)
            print('[INFO] Iteration ' + str(step + 1) + '. Base window length: ' + str(step_param.nfft))
            ## Starting spectral analysis
            # Segmentation/tapering
            station.apply_taper(step_param)
            # Robust regression
            station.recover_Z(step, project.parameters)
            if step < project.parameters.nb_reductions:
                print('[INFO] Iteration ' + str(step + 1) + ' done. Reducing base window length to: ' + str(int(step_param.nfft / (project.parameters.length_reduction))))
            print('[INFO] Response function computation done.')
        print("[INFO_TIME] Elapsed time for robust regression: " + str(time.time() - regression_time) + ' seconds.')
        print('[INFO] Now writing output file.')
        station.write_tensor(project)
        station.tensor.plot_tensor()
        if len(station.output) > 1: # 
            station.tensor.plot_tensor_res_phase()
