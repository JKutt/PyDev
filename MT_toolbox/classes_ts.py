import numpy as np
import pylab as plt
import spectrum_lib
import sys
import glob
import io_lib
from classes_z import Z
import multiprocessing


class project:

    def __init__(self,
                station=None,
                parameters=None,
                client=None,
                processor=None,
                version=None,
                folder=None,
                date_creation=None,
                type=None,
                output_level=None,
                log_file=None):

        self.station = station
        self.parameters = parameters
        self.client = client
        self.processor = processor
        self.version = version
        self.folder = folder
        self.date_creation = date_creation
        self.type = type # Z, T or Full
        self.output_level = output_level
        self.log_file = log_file

    def read_parameters_file(self):
        f_in = open('parameters.inp', 'r')
        lines = f_in.readlines()
        for line in lines:
            if 'station' in line:
                self.station = line.split('station=')[1].split('>')[0]
            if 'client' in line:
                self.client = line.split('client=')[1].split('>')[0]
            if 'processor' in line:
                self.processor = line.split('processor=')[1].split('>')[0]
            if 'folder' in line:
                self.folder = line.split('folder=')[1].split('>')[0]
            if 'type' in line:
                self.type = line.split('type=')[1].split('>')[0]
            if 'output_level' in line:
                self.output_level = int(line.split('output_level=')[1].split('>')[0])
        if self.station == None or\
           self.client == None or\
           self.processor == None or\
           self.folder == None or\
           self.type == None or\
           self.output_level == None:
            print('[ERROR] Missing fields in parameters file in the project section.')
            sys.exit()
        else:
            if self.output_level > 0:
                print('[INFO] Project section is loaded')
        self.parameters = parameters()
        for line in lines:
            if 'tbw' in line:
                self.parameters.tbw = int(line.split('tbw=')[1].split('>')[0])
            if 'nfft' in line:
                self.parameters.nfft = int(line.split('nfft=')[1].split('>')[0])
            if 'overlap' in line:
                self.parameters.overlap = int(line.split('overlap=')[1].split('>')[0])
            if 'taper' in line:
                self.parameters.taper = line.split('taper=')[1].split('>')[0]
            if 'index_first_frequency' in line:
                self.parameters.index_first_frequency = int(line.split('index_first_frequency=')[1].split('>')[0])
            if 'frequency_increment' in line:
                self.parameters.frequency_increment = int(line.split('frequency_increment=')[1].split('>')[0])
            if 'nb_increment' in line:
                self.parameters.nb_increment = int(line.split('nb_increment=')[1].split('>')[0])
            if 'length_reduction' in line:
                self.parameters.length_reduction = int(line.split('length_reduction=')[1].split('>')[0])
            if 'nb_reductions' in line:
                self.parameters.nb_reductions = int(line.split('nb_reductions=')[1].split('>')[0])
            if 'output=' in line:
                self.parameters.output = line.split('output=')[1].split('>')[0]
            if 'error=' in line:
                self.parameters.error = line.split('error=')[1].split('>')[0]

        if self.output_level > 0:
            print('[INFO] Parameters section is loaded')
        f_in.close()

    def check_project(self):
        var_ex = 0
        tmp = glob.glob(self.folder)
        if not len(tmp):
            print('[ERROR] Check project folder. Seems like it doesn\'t exists')
            var_ex = 1

        if not self.type in ['Z', 'T', 'Full']:
            print('[ERROR] Check project type.')
            var_ex = 1

        if not self.parameters.taper in ['slepian', 'hamming', 'window_tukey', 'null']:
            print('[ERROR] Check taper type.')
            var_ex = 1

        if var_ex:
            print('Read error messages and correct.')
            sys.exit()
        else:
            print('[INFO] Everything seems good to start reading files. Starting to write log file.')
            self.log_file = self.folder + '/processing.log'


class station_mt:

    def __init__(self,
                input=None,
                output=None,
                ref=None,
                tensor=None,
                output_level=None,
                sample_freq=None,
                log_file=None):

        self.input = []
        self.output = []
        self.ref = []
        self.tensor = Z()
        self.output_level = output_level
        if self.output_level == None:
            self.output_level = 0
        self.sample_freq = sample_freq
        self.log_file = log_file

    def read_parameters_file_and_load_data(self):
        f_in = open('parameters.inp', 'r')
        lines = f_in.readlines()
        for line in lines:
            if '<ex=' in line:
                tmp = line.split('<ex=')[1].split('>')[0]
                self.output.append(time_serie(file_name=tmp,
                                              sample_freq=2048,
                                              output_level=self.output_level,
                                              log_file=self.log_file))
            if '<ey' in line:
                tmp = line.split('<ey=')[1].split('>')[0]
                self.output.append(time_serie(file_name=tmp,
                                              sample_freq=2048,
                                              output_level=self.output_level,
                                              log_file=self.log_file))
            if '<hx' in line:
                tmp = line.split('<hx=')[1].split('>')[0]
                self.input.append(time_serie(file_name=tmp,
                                             sample_freq=2048,
                                             output_level=self.output_level,
                                             log_file=self.log_file))
            if '<hy' in line:
                tmp = line.split('<hy=')[1].split('>')[0]
                self.input.append(time_serie(file_name=tmp,
                                             sample_freq=2048,
                                             output_level=self.output_level,
                                             log_file=self.log_file))
            if '<rx' in line:
                tmp = line.split('<rx=')[1].split('>')[0]
                self.ref.append(time_serie(file_name=tmp,
                                           sample_freq=2048,
                                           output_level=self.output_level,
                                           log_file=self.log_file))
            if '<ry' in line:
                tmp = line.split('<ry=')[1].split('>')[0]
                self.ref.append(time_serie(file_name=tmp,
                                           sample_freq=2048,
                                           output_level=self.output_level,
                                           log_file=self.log_file))

        if self.output_level > 0:
            print('[INFO] Data section is loaded')
        f_in.close()

    def load_frequencies(self, params):
        self.tensor.frequencies = spectrum_lib.get_frequencies(self.sample_freq, params)
        self.tensor.Z = np.asarray([[[np.complex(0, 0), np.complex(0, 0)]
                                     for f in range(len(self.tensor.frequencies))]
                                     for channel in range(len(self.output))])

        self.tensor.error = np.asarray([[[np.complex(0, 0), np.complex(0, 0)]
                                          for f in range(len(self.tensor.frequencies))]
                                          for channel in range(len(self.output))])
        if self.output_level > 1:
            print('Analyzed frequencies:')
            for i, f in zip(range(len(self.tensor.frequencies)), self.tensor.frequencies):
                print('\t ' + str(i) + '. ' + str(f))


    def check_station(self, type):

        sample_freqs = []
        for i in range(len(self.output)):
            sample_freqs.append(self.output[i].sample_freq)
        for i in range(len(self.input)):
            sample_freqs.append(self.input[i].sample_freq)
        for i in range(len(self.ref)):
            sample_freqs.append(self.ref[i].sample_freq)

        if type == 'Z':
            if not len(self.output) == 2 or not len(self.input) == 2:
                print('[ERROR] Check type/input time series.')
                sys.exit()

        if type == 'T':
            if not len(self.output) == 1 or not len(self.input) == 2:
                print('[ERROR] Check type/input time series.')
                sys.exit()

        if type == 'Full':
            if not len(self.output) == 3 or not len(self.input) == 2:
                print('[ERROR] Check type/input time series.')
                sys.exit()

        sizes = []
        for i in range(len(self.output)):
            sizes.append(len(self.output[i].data))
        for i in range(len(self.input)):
            sizes.append(len(self.input[i].data))
        for i in range(len(self.ref)):
            sizes.append(len(self.ref[i].data))

        if len(np.unique(sample_freqs)) > 1:
            print('[ERROR] Sample frequencies are not equal. I am not capable of fixing this by myself for now.')
            sys.exit()
        else:
            self.sample_freq = sample_freqs[0]
            if self.output_level > 0:
                print('[INFO] Time series have equal sample frequencies.')

        if np.min(sizes) == 0:
            print('[ERROR] One of the time series is empty.')
            sys.exit()
        else:
            if np.max(sizes) > 0:
                if len(np.unique(sizes)) > 1:
                    print('[WARNING] Time series length are not equal. Code will only use the smallest length. Changing size now.')
                    for i in range(len(self.output)):
                        self.output[i].data = self.output[i].data[:np.min(sizes)]
                        self.output[i].length = int(np.min(sizes))
                    for i in range(len(self.input)):
                        self.input[i].data = self.input[i].data[:np.min(sizes)]
                        self.input[i].length = int(np.min(sizes))
                    for i in range(len(self.ref)):
                        self.ref[i].data = self.ref[i].data[:np.min(sizes)]
                        self.ref[i].length = int(np.min(sizes))
                    print('[WARNING] New size is now: ' + str(np.min(sizes)))
                else:
                    if self.output_level > 0:
                        print('[INFO] Time series have the same length.')

        print('[INFO] Everything seems good to start spectral analysis.')

    def plot(self):

        fig, (ax0, ax1, ax2, ax3, ax4, ax5) = plt.subplots(6, 1, sharex=True)
        ax0.plot(self.output[0].data)
        ax1.plot(self.output[1].data)
        ax2.plot(self.input[0].data)
        ax3.plot(self.input[1].data)
        ax4.plot(self.ref[0].data)
        ax5.plot(self.ref[1].data)
        plt.show()


    def apply_taper(self, param):
        print("here")
        for i in range(len(self.output)):
            self.output[i].segment_data(param)
            self.output[i].taper_segments(param)
            self.output[i].get_fft_segments()

        for i in range(len(self.input)):
            self.input[i].segment_data(param)
            self.input[i].taper_segments(param)
            self.input[i].get_fft_segments()

        for i in range(len(self.ref)):
            self.ref[i].segment_data(param)
            self.ref[i].taper_segments(param)
            self.ref[i].get_fft_segments()


    def recover_Z(self, step, param):

        ### PARALLEL
        # Using a pool of worker to work on each frequency
        max_number_processes = multiprocessing.cpu_count()
        #pool = multiprocessing.Pool(max_number_processes - 20)
        pool = multiprocessing.Pool(2)
        r = []
        print("\tStarting workers pool: " + str(2) + " workers working on robust regressions.")
        for ind in range(param.nb_increment):
            #res = pool.map(datlib.read_node_library, )
            index = param.index_first_frequency + ind * param.frequency_increment
            for channel in range(len(self.output)):
                r.append(pool.apply_async(spectrum_lib.robust_regression,
                                            (self.output[channel].get_fourier_index(index),
                                            [self.input[i].get_fourier_index(index) for i in range(len(self.input))],
                                            [self.ref[i].get_fourier_index(index) for i in range(len(self.ref))],
                                            self.output_level)))
        pool.close()
        pool.join()

        # Getting results from multiprocessing
        for ind in range(param.nb_increment):
            for channel in range(len(self.output)):
                tmp = r[ind * len(self.output) + channel].get()
                self.tensor.Z[channel, step * param.nb_increment + ind, :] = tmp[0]
                self.tensor.error[channel, step * param.nb_increment + ind, :] = tmp[1]

        """
        ### SERIAL
        for ind in range(param.nb_increment):
            print('\t[INFO] Frequency increment ' + str(ind + 1) + '/' + str(param.nb_increment))
            index = param.index_first_frequency + ind * param.frequency_increment
            for channel in range(len(self.output)):
                #self.tensor.Z[channel, step * param.nb_increment + ind, :], \
                #self.tensor.error[channel, step * param.nb_increment + ind, :] =
                self.tensor.Z[channel, step * param.nb_increment + ind, :], \
                self.tensor.error[channel, step * param.nb_increment + ind, :] = spectrum_lib.robust_regression(
                                                self.output[channel].get_fourier_index(index),
                                                [self.input[i].get_fourier_index(index) for i in range(len(self.input))],
                                                [self.ref[i].get_fourier_index(index) for i in range(len(self.ref))],
                                                self.output_level
                                                                                                            )
        """
    def rotate_station(self, angle):
        print('TODO')



    def write_Z(self, project):
        print('TODO')


class parameters:

    def __init__(self,
                output=None,
                tbw=None,
                nfft=None,
                overlap=None,
                taper=None,
                index_first_frequency=None,
                frequency_increment=None,
                nb_increment=None,
                length_reduction=None,
                nb_reductions=None,
                output_level=None,
                error=None):
        self.output = output
        self.nfft = nfft
        self.tbw = tbw
        self.overlap = overlap
        self.taper = taper
        self.index_first_frequency = index_first_frequency
        self.frequency_increment = frequency_increment
        self.nb_increment = nb_increment
        self.length_reduction = length_reduction
        self.nb_reductions = nb_reductions
        self.error = error


class step_parameters:
    def __init__(self,
                tbw=None,
                nfft=None,
                overlap=None,
                taper=None,
                index_first_frequency=None,
                frequency_increment=None,
                nb_increment=None,
                error=None):

        self.nfft = nfft
        self.tbw = tbw
        self.overlap = overlap
        self.taper = taper
        self.index_first_frequency = index_first_frequency
        self.frequency_increment = frequency_increment
        self.nb_increment = nb_increment
        self.error = error

    def update_param(self, parameters, iteration):
        self.nfft = int(parameters.nfft / (parameters.length_reduction ** iteration))
        self.tbw = parameters.tbw
        self.overlap = parameters.overlap
        self.taper = parameters.taper
        self.index_first_frequency = parameters.index_first_frequency
        self.frequency_increment = parameters.frequency_increment
        self.nb_increment = parameters.nb_increment
        self.error = parameters.error


class time_serie:

    def __init__(self,
                file_name=None,
                date_start=None,
                sample_freq=None,
                length=None,
                easting=None,
                northing=None,
                elevation=None,
                type=None,
                orientation=None,
                data=None,
                data_fft=None,
                segments=None,
                segments_fft=None,
                output_level=None,
                log_file=None):
        self.file_name = file_name
        self.date_start = date_start
        self.sample_freq = sample_freq
        self.length = length
        self.easting = easting
        self.northing = northing
        self.elevation = elevation
        self.type = type
        self.orientation = orientation
        self.data = data
        self.data_fft = data_fft
        self.segments = []
        self.segments_fft = []
        self.output_level = output_level
        self.log_file = log_file
        self.check_path()
        self.load_data()

    def check_path(self):
        tmp = glob.glob(self.file_name)
        if not len(tmp):
            print('[ERROR] Path: ' + self.file_name + ' doesn\'t exist. Stopping.')
            sys.exit()
        else:
            if self.output_level > 0:
                print('[INFO] Path: ' + self.file_name + ' exists. Loading time series.')

    def load_data(self):
        try:
            self.data=io_lib.read_binary_file(self.file_name)
            #self.data=io_lib.read_ascii_file(self.file_name)
        except:
            print('[ERROR] Something went wrong while loading the file ' + self.file_name)
            sys.exit()

    def segment_data(self, parameters):
        self.segments = []
        n_overlap = int(np.floor(parameters.nfft - parameters.nfft * parameters.overlap / 100))
        n_windows = int(np.floor(len(self.data) / n_overlap)) - 1

        for ind in range(n_windows):
            self.segments.append(self.data[n_overlap * ind:n_overlap * ind + parameters.nfft])

    def taper_segments(self, parameters):
        if len(self.segments):
            taper = spectrum_lib.get_taper(parameters.nfft, parameters.taper, parameters.tbw)
            for i in range(len(self.segments)):
                spectrum_lib.apply_taper(self.segments[i], taper)
        else:
            print('[WARNING] Segmentation hasn\'t been done yet. Use self.segment_data(self, parameters).')

    def get_fft_segments(self):
        self.segments_fft = []
        for i in range(len(self.segments)):
            self.segments_fft.append(np.fft.fft(self.segments[i]))

    def get_fourier_index(self, index):
        array_spectral = np.asarray([np.complex(0, 0) for i in range(len(self.segments))])
        for i, segment in zip(range(len(self.segments)), self.segments_fft):
            array_spectral[i] = segment[index]

        return array_spectral
