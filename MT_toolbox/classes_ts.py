import numpy as np
import pylab as plt
import spectrum_lib
import wavelet_lib as wav_lib
import sys
import glob
import io_lib
from classes_z import Z
import multiprocessing


class project:

    def __init__(self,
                patch=None,
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

        self.patch = patch
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
            if 'patch' in line:
                self.patch = line.split('patch=')[1].split('>')[0]
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
        if self.patch == None or \
           self.station == None or\
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
            if 'local_coherence=' in line:
                self.parameters.local_coherence = float(line.split('local_coherence=')[1].split('>')[0])
            if 'remote_coherence=' in line:
                self.parameters.remote_coherence = float(line.split('remote_coherence=')[1].split('>')[0])
            if 'wavelets=' in line:
                self.parameters.wavelets = line.split('wavelets=')[1].split('>')[0]

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

        if self.parameters.remote_coherence == None:
            print('[INFO] No coherence thresholding between local and remote has been setup in the parameters files. Will not use.')

        if self.parameters.local_coherence == None:
            print('[INFO] No coherence thresholding between local output and local input has been setup in the parameters files. Will not use.')


class station_mt:

    def __init__(self,
                input=None,
                output=None,
                ref=None,
                tensor=None,
                output_level=None,
                sample_freq=None,
                log_file=None,
                location=None,
                coherence=None):

        self.input = []
        self.output = []
        self.ref = []
        self.tensor = Z()
        self.output_level = output_level
        if self.output_level == None:
            self.output_level = 0
        self.sample_freq = sample_freq
        self.log_file = log_file
        self.location = location
        self.coherence = coherence

    def read_parameters_file_and_load_data(self):
        f_in = open('parameters.inp', 'r')
        lines = f_in.readlines()
        for line in lines:
            if '<ex=' in line:
                tmp = line.split('<ex=')[1].split('>')[0]
                self.output.append(time_serie(file_name=tmp,
                                              sample_freq=2048,
                                              output_level=self.output_level,
                                              log_file=self.log_file,
                                              type='Ex'))
            if '<ey' in line:
                tmp = line.split('<ey=')[1].split('>')[0]
                self.output.append(time_serie(file_name=tmp,
                                              sample_freq=2048,
                                              output_level=self.output_level,
                                              log_file=self.log_file,
                                              type='Ey'))

            if '<hz' in line:
                tmp = line.split('<hz=')[1].split('>')[0]
                self.output.append(time_serie(file_name=tmp,
                                              sample_freq=2048,
                                              output_level=self.output_level,
                                              log_file=self.log_file,
                                              type='Hz'))
            if '<hx' in line:
                tmp = line.split('<hx=')[1].split('>')[0]
                self.input.append(time_serie(file_name=tmp,
                                             sample_freq=2048,
                                             output_level=self.output_level,
                                             log_file=self.log_file,
                                             type='Hx'))
            if '<hy' in line:
                tmp = line.split('<hy=')[1].split('>')[0]
                self.input.append(time_serie(file_name=tmp,
                                             sample_freq=2048,
                                             output_level=self.output_level,
                                             log_file=self.log_file,
                                             type='Hy'))
            if '<rx' in line:
                tmp = line.split('<rx=')[1].split('>')[0]
                self.ref.append(time_serie(file_name=tmp,
                                           sample_freq=2048,
                                           output_level=self.output_level,
                                           log_file=self.log_file,
                                           type='Rx'))
            if '<ry' in line:
                tmp = line.split('<ry=')[1].split('>')[0]
                self.ref.append(time_serie(file_name=tmp,
                                           sample_freq=2048,
                                           output_level=self.output_level,
                                           log_file=self.log_file,
                                           type='Ry'))

        if self.output_level > 0:
            print('[INFO] Data section is loaded')
        f_in.close()

    def load_frequencies(self, params):
        self.tensor.frequencies = spectrum_lib.get_frequencies(self.sample_freq, params)
        self.tensor.Z = np.asarray([[[np.complex(0, 0), np.complex(0, 0)]
                                     for f in range(len(self.tensor.frequencies))]
                                     for channel in range(len(self.output))])

        self.tensor.error = np.asarray([[[0., 0.]
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
        self.location = [self.output[0].northing, self.output[0].easting,
                         self.output[0].elevation, self.output[0].zone]

    def get_wavelet_transforms(self):
        print('[WARNING] Hugo: Needs to parametrize input frequencies. Will do it later.')
        dyadic_frequencies = wav_lib.dyadic_frequencies(2, 64, 10)
        for i in range(len(self.output)):
            self.output[i].wavelets = wav_lib.wavelet_analysis(self.output[i].data, self.output[i].sample_freq, dyadic_frequencies)
        for i in range(len(self.input)):
            self.input[i].wavelets = wav_lib.wavelet_analysis(self.input[i].data, self.input[i].sample_freq, dyadic_frequencies)
            fig, ax1 = plt.subplots(1, 1)
            ax1.imshow(self.input[i].wavelets.transpose(), aspect='auto', cmap='jet')
            plt.show()
        for i in range(len(self.ref)):
            self.ref[i].wavelets = wav_lib.wavelet_analysis(self.ref[i].data, self.ref[i].sample_freq, dyadic_frequencies)

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


    def get_coherence_elec_mag(self, index):

        self.coherence = [[], []]
        for j in range(len(self.input)):
            for i in range(len(self.output[0].segments)):
                coherence = spectrum_lib.get_coherence_window(self.output[j].segments[i],
                                            self.input[j].segments[i],
                                            self.input[j].sample_freq)
                self.coherence[j].append(coherence)


        array_coherence = np.asarray([[0. for i in range(len(self.output[0].segments))] for j in range(2)])
        for i, segment in zip(range(len(self.output[0].segments)), self.coherence[0]):
            array_coherence[0][i] = self.coherence[0][i][index]
            array_coherence[1][i] = self.coherence[1][i][index]

        return array_coherence

    def get_coherence_local_remote(self, index):

        self.coherence = [[], []]
        for j in range(len(self.input)):
            for i in range(len(self.ref[0].segments)):
                coherence = spectrum_lib.get_coherence_window(self.input[j].segments[i],
                                            self.ref[j].segments[i],
                                            self.input[j].sample_freq)
                self.coherence[j].append(coherence)


        array_coherence = np.asarray([[0. for i in range(len(self.ref[0].segments))] for j in range(2)])
        for i, segment in zip(range(len(self.ref[0].segments)), self.coherence[0]):
            array_coherence[0][i] = self.coherence[0][i][index]
            array_coherence[1][i] = self.coherence[1][i][index]

        return array_coherence

    def apply_coherence(self, output_in, input_in, ref_in, array_coherence, criteria):

        output = []
        input = [[] for i in range(len(input_in))]
        ref = [[] for i in range(len(ref_in))]
        for i in range(len(output_in)):
            #if not np.isnan(array_coherence[0][i]) and not np.isnan(array_coherence[1][i]):
            if array_coherence[0][i] > criteria or array_coherence[1][i] > criteria:
                output.append(output_in[i])
                for j in range(len(input_in)):
                    input[j].append(input_in[j][i])
                for j in range(len(ref_in)):
                    ref[j].append(ref_in[j][i])

        return output, input, ref

    def recover_Z(self, step, param):

        ### PARALLEL
        # Using a pool of worker to work on each frequency

        max_number_processes = multiprocessing.cpu_count()
        #pool = multiprocessing.Pool(max_number_processes - 20)
        pool = multiprocessing.Pool(2)
        r = []
        
        print("\tStarting workers pool: " + str(2) + " workers working on robust regressions.")
        for ind in range(param.nb_increment):
            index = param.index_first_frequency + ind * param.frequency_increment
            print("\tApplying coherence thresholding.")
            if not param.remote_coherence == None:
                array_coherence_remote = self.get_coherence_local_remote(index)
            if not param.local_coherence == None:
                array_coherence_local = self.get_coherence_elec_mag(index)
            for channel in range(len(self.output)):
                if not param.remote_coherence == None:
                    output_reg, input_reg, ref_reg = self.apply_coherence(
                                            self.output[channel].get_fourier_index(index),
                                            [self.input[i].get_fourier_index(index) for i in range(len(self.input))],
                                            [self.ref[i].get_fourier_index(index) for i in range(len(self.ref))],
                                            array_coherence_remote, param.remote_coherence)
                else:
                    output_reg = self.output[channel].get_fourier_index(index)
                    input_reg = [self.input[i].get_fourier_index(index) for i in range(len(self.input))]
                    ref_reg = [self.ref[i].get_fourier_index(index) for i in range(len(self.ref))]
                if not param.local_coherence == None:
                    output_reg, input_reg, ref_reg = self.apply_coherence(
                                            output_reg, input_reg, ref_reg, array_coherence_local, param.local_coherence)

                r.append(pool.apply_async(spectrum_lib.robust_regression,
                                            (output_reg, input_reg, ref_reg, self.output_level)))
                #r.append(pool.apply_async(spectrum_lib.robust_regression,
                #                            (self.output[channel].get_fourier_index(index),
                #                            [self.input[i].get_fourier_index(index) for i in range(len(self.input))],
                #                            [self.ref[i].get_fourier_index(index) for i in range(len(self.ref))],
                #                            self.output_level)))
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



    def write_tensor(self, project):
        fOut = open(project.folder + '/' + project.station + '.DAT', 'w')
        fOut.write('#Patch: ' + project.patch + '\n')
        fOut.write('#Station: ' + project.station + '\n')
        if not None in self.location:
            fOut.write('#Location (UTM): ' + '{:10f}'.format(self.location[0]) + ';'
                    '{:10f}'.format(self.location[1]) + ';'
                    '{:10f}'.format(self.location[2]) + ';Zone:' +
                    self.location[3] + '\n')
        else:
            fOut.write('#Location (UTM): NA;NA;NA;Zone:NA' + '\n')
        fOut.write('#Client: ' + project.client + '\n')
        fOut.write('#Processor: ' + project.processor + '\n')
        fOut.write('#Time series:\n')
        for ts in self.output:
            fOut.write('#' + ts.type + ': ' + ts.file_name + '\n')
        for ts in self.input:
            fOut.write('#' + ts.type + ': ' + ts.file_name + '\n')
        for ts in self.ref:
            fOut.write('#' + ts.type + ': ' + ts.file_name + '\n')

        self.tensor.get_res_phase_error()
        # Writing real and imaginary parts
        if project.type == 'Z':
            fOut.write('#Tensor type: Impedance (Real and Imaginary)\n')
            fOut.write('#Units: Frequency [Hz], Z [mV/km/nT]\n')
            columns = ['#Frequency', 'R Zxx',
                        'R Zxy', 'R Zyx',
                        'R Zyy', 'I Zxx',
                        'I Zxy', 'R Zyxr',
                        'I Zyy', 'STD Zxx',
                        'STD Zxy', 'STD Zyx',
                        'STD Zyy']
            for item in columns:
                fOut.write('{:>15}'.format(item))
            fOut.write('\n')
        if project.type == 'T':
            fOut.write('Tensor type: Tipper\n')
            fOut.write('#Units: T []\n')
            columns = ['#Frequency', 'R Tx',
                        'I Tx', 'R Ty',
                        'I Ty', 'STD Tx',
                        'STD Ty']
            for item in columns:
                fOut.write('{:>15}'.format(item))
            fOut.write('\n')
        if project.type == 'Full':
            fOut.write('Tensor type: Full tensor (Impedance + Tipper)\n')
            fOut.write('#Units: Frequency [Hz], Z [mV/km/nT], T []\n')
            columns = ['#Frequency', 'R Zxx',
                        'R Zxy', 'R Zyx',
                        'R Zyy', 'R Tx',
                        'R Ty', 'I Zxx',
                        'I Zxy', 'R Zyxr',
                        'I Zyy', 'I Tx',
                        'I Ty', 'STD Zxx',
                        'STD Zxy', 'STD Zyx',
                        'STD Zyy', 'STD Tx',
                        'STD Ty']
            for item in columns:
                fOut.write('{:>15}'.format(item))
            fOut.write('\n')
        for (freq, ind) in zip(self.tensor.frequencies, range(len(self.tensor.frequencies))):
            fOut.write('{:15f}'.format(freq))
            for channel in range(len(self.tensor.Z)):
                fOut.write('{:15f}'.format(np.real(self.tensor.Z[channel, ind, 0])))
                fOut.write('{:15f}'.format(np.real(self.tensor.Z[channel, ind, 1])))
            for channel in range(len(self.tensor.Z)):
                fOut.write('{:15f}'.format(np.imag(self.tensor.Z[channel, ind, 0])))
                fOut.write('{:15f}'.format(np.imag(self.tensor.Z[channel, ind, 1])))

            for channel in range(len(self.tensor.error)):
                fOut.write('{:15f}'.format((self.tensor.error[channel, ind, 0])))
                fOut.write('{:15f}'.format((self.tensor.error[channel, ind, 1])))
            fOut.write('\n')

        # Writing rho/phase parts
        if project.type == 'Z' or project.type == 'Full':
            fOut.write('#Tensor type: Impedance (Resistivity and Phase)\n')
            fOut.write('#Units: Frequency [Hz], Z [mV/km/nT]\n')
            columns = ['#Frequency', 'Rho_xx',
                        'Rho_xy', 'Rho_yx',
                        'Rho_yy', 'Phase_xx',
                        'Phase_xy', 'Phase_yx',
                        'Phase_yy', 'STD Rho_xx',
                        'STD Rho_xy', 'STD Rho_yx',
                        'STD Rho_yy', 'STD Phase_xx',
                        'STD Phase_xy', 'STD Phase_yx',
                        'STD Phase_yy']
            for item in columns:
                fOut.write('{:>15}'.format(item))
            fOut.write('\n')

            for (freq, ind) in zip(self.tensor.frequencies, range(len(self.tensor.frequencies))):
                fOut.write('{:15f}'.format(freq))
                for channel in range(2):
                    fOut.write('{:15f}'.format(np.abs(self.tensor.Z[channel, ind, 0]) ** 2 / 5. / freq))
                    fOut.write('{:15f}'.format(np.abs(self.tensor.Z[channel, ind, 1]) ** 2 / 5. / freq))
                for channel in range(2):
                    fOut.write('{:15f}'.format(np.angle(self.tensor.Z[channel, ind, 0])))
                    fOut.write('{:15f}'.format(np.angle(self.tensor.Z[channel, ind, 1])))

                for channel in range(2):
                    fOut.write('{:15f}'.format(self.tensor.res_error[channel, ind, 0]))
                    fOut.write('{:15f}'.format(self.tensor.res_error[channel, ind, 1]))
                for channel in range(2):
                    fOut.write('{:15f}'.format(self.tensor.phase_error[channel, ind, 0]))
                    fOut.write('{:15f}'.format(self.tensor.phase_error[channel, ind, 1]))
                fOut.write('\n')
        fOut.close()

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
                error=None,
                local_coherence=None,
                remote_coherence=None,
                wavelets=None):
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
        self.local_coherence = local_coherence
        self.remote_coherence = remote_coherence
        self.wavelets = wavelets


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
                zone=None,
                type=None,
                orientation=None,
                data=None,
                data_fft=None,
                segments=None,
                segments_fft=None,
                segments_coherence=None,
                wavelets=None,
                output_level=None,
                log_file=None):
        self.file_name = file_name
        self.date_start = date_start
        self.sample_freq = sample_freq
        self.length = length
        self.easting = easting
        self.northing = northing
        self.elevation = elevation
        self.zone = zone
        self.type = type
        self.orientation = orientation
        self.data = data
        self.data_fft = data_fft
        self.segments = []
        self.segments_fft = []
        self.segments_coherence = []
        self.wavelets = wavelets
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
