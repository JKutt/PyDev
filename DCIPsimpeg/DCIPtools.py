#################################################
# Imports

import numpy as np
from scipy import fftpack
# from scipy import sparse
import matplotlib.pyplot as plt
from scipy import signal
from scipy.special import factorial
from SimPEG.EM.Static import DC
import properties
from numpy import genfromtxt
##################################################
# DCIPtools Version 1.0  April 2019


class baseKernel(object):

    kernel = properties.Array(
        "array containing filter kernal",
        default=np.zeros(1)
    )

    num_ensembles = properties.Integer(
        "number of ensembles",
        default=0)

    reject_std = properties.Float(
        "Tolerance for rejecting an ensemble or stat reject",
        default=1.,
        required=True)

    def __init__(self, filtershape=None, **kwargs):
        if filtershape is not None:
            if isinstance(filtershape, np.ndarray):
                self.filtershape = filtershape
                self.createFilterKernal()

    def createFilterKernal(self):
        """
        creates the filter kernal

        """
        if len(self.filtershape) > 1:
            # create the filter kernal
            tkernal = np.ones((3, 1))
            tkernal[1] = -2.0                    # 3 point kernal
            bsgn = np.ones((1, self.filtershape.size))
            bsgn[0, 1::2] = bsgn[0, 1::2] * -1   # alternate pol rem lin drift
            bwtd = np.matmul(tkernal,
                             bsgn * self.filtershape)  # filter weights
            tmp1 = np.arange(1, 4)
            tmp1 = np.reshape(tmp1, (3, 1))
            tmp2 = np.ones((1, self.filtershape.size))
            tmp3 = np.ones((3, 1))
            tmp4 = np.arange(self.filtershape.size)
            tmp4 = np.reshape(tmp4, (1, self.filtershape.size))
            knew = (np.matmul(tmp1, tmp2) +
                    np.matmul(tmp3, (tmp4 * (self.filtershape.size + 3))))
            print(knew)
            btmp = np.zeros((self.filtershape.size + 2,
                             self.filtershape.size))  # create zero matrix
            shape_knew = knew.shape
            num_elements_kn = shape_knew[0] * shape_knew[1]
            knew = np.reshape(knew, num_elements_kn, order='F')
            shape = btmp.shape
            num_elements = shape[0] * shape[1]
            btmp = np.reshape(btmp, num_elements, order='F')
            shape_bwtd = bwtd.shape
            num_elements_b = shape_bwtd[0] * shape_bwtd[1]
            bwtd = np.reshape(bwtd, num_elements_b, order='F')
            for idx in range(knew.size):
                btmp[int(knew[idx]) - 1] = bwtd[idx]  # fill diag w/ weights
            btmp = np.reshape(btmp, shape, order='F')
            tHK = np.sum(btmp, 1)
            norm_tHK = np.sum(np.abs(tHK))
            tHK = tHK / norm_tHK
            tHK = np.reshape(tHK, (tHK.size, 1))
        else:
            raise Exception('filter size must be greater than 3!')
            tHK = np.zeros(1)
        self.kernel = tHK[:self.filtershape.size]         # assign weighted kernal

    def __mul__(self, val):
        """
        stack results

        """
        if isinstance(val, np.ndarray):
            return self._transform(val)
        else:
            raise Exception("Input must be an numpy array")

    @property
    def sizeOfFilter(self):
        """
            :rtypr int
            :return: number of points in filter kernal
        """
        return self.filtershape.size

    @property
    def getFilterKernal(self):
        """
           returns the filter kernal
        """
        return self.kernel


class decayKernal(baseKernel):
    """
       Decay kernal for calculating a decay
       from a stack
    """
    # output_type = properties.StringChoice(
    #     "tells kernal the type of output (decay, std)",
    #     default="Vs",
    #     choices=["std", "Vs"]
    # )

    # window_starts = properties.Array(
    #     "array containing window start times",
    #     default=np.zeros(20)
    # )

    # window_ends = properties.Array(
    #     "array containing window end times",
    #     default=np.zeros(1)
    # )

    # window_weight = properties.Float(
    #     "attenuation",
    #     default=401,
    #     required=True
    # )

    # window_overlap = properties.Float(
    #     "overlap in windows",
    #     default=0.45,
    #     required=True
    # )

    def __init__(self, num_windows=None,
                 window_starts=None,
                 window_ends=None,
                 window_weight=None,
                 window_overlap=None,
                 output_type=None, **kwargs):
        baseKernel.__init__(self, None, **kwargs)

        if num_windows is not None:
            self.num_window = int(num_windows)
            self.window_starts = window_starts
            self.window_ends = window_ends
            self.window_weight = window_weight
            self.window_overlap = window_overlap
            self.output_type = output_type
        else:
            raise Exception("need number of windows: num_windows")
        # self.window_overlap = self.window_overlap / self.divsor

    def _transform(self, stack):
        """
        takes in stack data and returns decay
        Input:
        stack = Half period stacked Voltage data

        """
        # calculate weighted window decay data =============
        if isinstance(stack, np.ndarray):
            print(self.output_type)
            starts = self.getWindowStarts()
            # print(starts.size)
            timebase = np.round(
                starts[self.num_window - 1] / 1000.)
            timebase = timebase * 1000
            time = np.arange(0, timebase, (timebase / stack.size))
            vsDecay = np.zeros((self.num_window))
            vs_std = np.zeros((self.num_window))
            # loop through and do every window
            for win in range(self.num_window):
                self.window_starts[win]
                # find how many samples in first window
                cntr = 0
                # get time span for tail of windows
                start_tmp = self.window_starts[win] - (
                    self.window_overlap *
                    (self.window_ends[win] - self.window_starts[win]))
                if win == (vsDecay.size - 1):
                    end_tmp = self.window_ends[win]
                else:
                    end_tmp = self.window_ends[win] + (
                        self.window_overlap *
                        (self.window_ends[win] - self.window_starts[win]))
                # print start_tmp, end_tmp
                for i in range(stack.size):
                    if time[i] >= start_tmp and time[i] <= end_tmp:
                        cntr += 1
                # create window wieghts
                indx1 = np.arange(0, self.window_weight)
                weights = 0.5 - (0.5 *
                                 np.cos((2 * np.pi * indx1) /
                                        (indx1.size - 1)))
                # create new weights
                Wg = np.zeros(cntr)
                start_Wg = (indx1.size / 2.0 - 1.0) - (cntr / 2.0) + 1
                for r in range(Wg.size):
                    Wg[r] = weights[int(start_Wg) + r]
                # print Wg
                # create vector storing weighted values
                Vs_window = np.zeros(cntr)
                Vs_window_ave = np.zeros(cntr)
                # get window times
                w_idx = np.zeros(cntr)
                # assign total time and time step
                count = 0
                for i in range(time.size):
                    if time[i] >= start_tmp and time[i] <= end_tmp:
                        w_idx[count] = time[i]
                        Vs_window[count] = stack[i] * -1 * Wg[count]
                        Vs_window_ave[count] = stack[i] * -1
                        count += 1
                sumWin = np.sum(Vs_window)      # sum the values of the window
                vs_std[win] = np.std(Vs_window_ave)  # standard deviation of window
                # print Vs_window
                vsDecay[win] = sumWin / cntr
        else:
            raise Exception("input must be a stack numpy array!")
            vsDecay = np.zeros(1)
            # end decay =======================================
        output = vsDecay
        # output = vs_std
        if self.output_type == 'std':
            output = vs_std

        return output

    def getWindowStarts(self):
        return self.window_starts

    def getWindowCenters(self):
        """
           returns window centers
        """
        window_centers = (np.asarray(self.window_ends) +
                          np.asarray(self.window_starts)) / 2.
        return window_centers

    def getWindowWidths(self):
        """
           returns window widths
        """
        window_widths = (np.asarray(self.window_ends) -
                         np.asarray(self.window_starts))
        return window_widths


class frequencyWindowingKernal(baseKernel):
    """
       Frequency windowing kernal for evaluating frequencies
       from a time-series recording
    """

    def __init__(self,
                 eval_frequencies=None,
                 window_overlap=None,
                 output_type=None, **kwargs):
        baseKernel.__init__(self, None, **kwargs)

        if eval_frequencies is not None:
            self.eval_frequencies = eval_frequencies
            self.num_window = self.eval_frequencies.size
            self.window_overlap = window_overlap
            self.output_type = output_type
        else:
            raise Exception("need eval. frequencies: eval_frequencies")
        # self.window_overlap = self.window_overlap / self.divsor

    def _transform(self, time_series):
        """
        takes in stack data and returns decay
        Input:
        stack = Half period stacked Voltage data

        """
        # calculate weighted window decay data =============
        if isinstance(time_series, np.ndarray):
            print(self.output_type)
            if time_series.shape.size > 1:         # check mxn - stack results
                print("stacking similar frequency windows")
            else:                                  # continue with windowing
                print("windowing a single time-series")
            # starts = self.getWindowStarts()
            # # create window wieghts
            # indx1 = np.arange(0, self.window_weight)
            # # slepian goes here
            # window = signal.slepian(51, width=0.3)
            # weights = 0.5 - (0.5 *
            #                  np.cos((2 * np.pi * indx1) /
            #                         (indx1.size - 1)))
            # # create new weights
            # Wg = np.zeros(cntr)
            # start_Wg = (indx1.size / 2.0 - 1.0) - (cntr / 2.0) + 1
            # for r in range(Wg.size):
            #     Wg[r] = weights[int(start_Wg) + r]
            # # print Wg
            # # create vector storing weighted values
            # Vs_window = np.zeros(cntr)
            # Vs_window_ave = np.zeros(cntr)
            # # get window times
            # w_idx = np.zeros(cntr)
            # # assign total time and time step
            # count = 0
            # for i in range(time.size):
            #     if time[i] >= start_tmp and time[i] <= end_tmp:
            #         w_idx[count] = time[i]
            #         Vs_window[count] = stack[i] * -1 * Wg[count]
            #         Vs_window_ave[count] = stack[i] * -1
            #         count += 1
            # sumWin = np.sum(Vs_window)      # sum the values of the window
            # vs_std[win] = np.std(Vs_window_ave)  # standard deviation of window
            # # print Vs_window
            # vsDecay[win] = sumWin / cntr
        else:
            raise Exception("input must be a stack numpy array!")
            vsDecay = np.zeros(1)
            # end decay =======================================
        output = vsDecay
        # output = vs_std
        if self.output_type == 'std':
            output = 0 #vs_std

        return output

    def getWindowStarts(self):
        return self.window_starts

    def getWindowWidths(self):
        """
           returns window widths
        """
        # window_widths = (np.asarray(self.window_ends) -
        #                  np.asarray(self.window_starts))
        return window_widths


class filterKernal(baseKernel):
    """
        Filter Kernal for stacking a time-series
        raw signal
    """

    def __init__(self, filtershape, **kwargs):
        baseKernel.__init__(self, filtershape, **kwargs)
        # self.createFilterKernal()            # create filter kernal

    def _transform(self, signal):
        """
           performs the stacking calculation
        """
        size_of_stack = int(signal.size / (self.kernel.size))
        print("kernal size: {0} & Vt size: {1}".format(self.kernel.size, signal.size))
        Ax = np.reshape(signal, (int(size_of_stack),
                        int(self.kernel.size)), order='F')
        shape_Ax = Ax.shape
        shape_tHK = self.kernel.shape
        if shape_Ax[1] == shape_tHK[0]:
            stack = np.matmul(Ax, self.kernel)  # create stack data
            return stack
        else:
            return 0


class ensembleKernal(baseKernel):
    """
    ensemble stacking kernel

    """

    def __init__(self, filtershape,
                 number_half_periods, **kwargs):
        filter_kernal = createHanningWindow(filtershape.size)   # mod for tappered overlap
        baseKernel.__init__(self, filter_kernal, **kwargs)          # create the base kernal
        self.number_half_periods = number_half_periods              # number of half t in signal
        self.kernal_ends = self.createFilterKernalEnds(filtershape.size)  # create end kernal

    def createFilterKernalEnds(self, ensemble_size):
        """
        creates the filter kernal

        """
        ends_window = createHanningWindow(ensemble_size - 2)
        print(ends_window.size)
        if ends_window.size > 1:
            # create the filter kernal
            tkernal = np.ones((3, 1))
            tkernal[1] = -2.0                    # 3 point kernal
            bsgn = np.ones((1, ends_window.size))
            bsgn[0, 1::2] = bsgn[0, 1::2] * -1   # alternate pol rem lin drift
            bwtd = np.matmul(tkernal,
                             bsgn * ends_window)  # filter weights
            tmp1 = np.arange(1, 4)
            tmp1 = np.reshape(tmp1, (3, 1))
            tmp2 = np.ones((1, ends_window.size))
            tmp3 = np.ones((3, 1))
            tmp4 = np.arange(ends_window.size)
            tmp4 = np.reshape(tmp4, (1, ends_window.size))
            knew = (np.matmul(tmp1, tmp2) +
                    np.matmul(tmp3, (tmp4 * (ends_window.size + 3))))
            btmp = np.zeros((ends_window.size + 2,
                             ends_window.size))  # create zero matrix
            shape_knew = knew.shape
            num_elements_kn = shape_knew[0] * shape_knew[1]
            knew = np.reshape(knew, num_elements_kn, order='F')
            shape = btmp.shape
            num_elements = shape[0] * shape[1]
            btmp = np.reshape(btmp, num_elements, order='F')
            shape_bwtd = bwtd.shape
            num_elements_b = shape_bwtd[0] * shape_bwtd[1]
            bwtd = np.reshape(bwtd, num_elements_b, order='F')
            for idx in range(knew.size):
                btmp[int(knew[idx]) - 1] = bwtd[idx]  # fill diag w/ weights
            btmp = np.reshape(btmp, shape, order='F')
            tHK = np.sum(btmp, 1)
            norm_tHK = np.sum(np.abs(tHK))
            tHK = tHK / norm_tHK
            tHK = np.reshape(tHK, (tHK.size, 1))
        else:
            raise Exception('filter size must be greater than 3!')
            tHK = np.zeros(1)
        return tHK

    def _transform(self, signal):
        """
           performs the stacking calculation using Ensembles
        """
        sub_signals = []
        sub_samples = []
        size_of_stack = int(signal.size / (self.number_half_periods))
        overlap = 2                                                        # hard code standard
        T_per_ensemble = (self.kernel.size - overlap)                 # desired half T
        new_signal_size = int((T_per_ensemble) * size_of_stack)            # size of the signal going in
        samples_per_ens = int((T_per_ensemble + overlap) * size_of_stack)
        number_of_ensembles = signal.size / (T_per_ensemble * size_of_stack)
        opt_number_of_samples = number_of_ensembles * T_per_ensemble * size_of_stack
        print("num of ensembles: {0} & opt num: {1} & size of org sig: {2}".format(number_of_ensembles, T_per_ensemble, self.kernal_ends.size))
        signal = signal[:opt_number_of_samples]
        ensembles = np.zeros((size_of_stack, number_of_ensembles))         # matrix holding all the stacks 
        for index in range(number_of_ensembles):
            if index == 0:
                T_overlap = T_per_ensemble                                # get end overlap
                end_index = T_overlap * size_of_stack
                trim_signal = signal[:end_index]
                sub_signals.append(trim_signal)
                sub_samples.append(np.arange(0, end_index))
                Ax = np.reshape(trim_signal, (int(size_of_stack),
                                int(self.kernal_ends.size)), order='F')
                shape_Ax = Ax.shape
                shape_tHK = self.kernal_ends.shape
                if shape_Ax[1] == shape_tHK[0]:
                    stack = np.matmul(Ax, self.kernal_ends)  # create stack data
                    ensembles[:, index] = stack.T
                else:
                    print("fail stack, wrong size")
            elif index == 1:
                T_overlap = T_per_ensemble + overlap
                start_index = size_of_stack * (T_per_ensemble - 2)
                end_index = start_index + (T_overlap * size_of_stack)
                trim_signal = signal[start_index:end_index]
                sub_signals.append(trim_signal)
                sub_samples.append(np.arange(start_index, end_index))
                print("num of ensembles: {0} & opt num: {1} & size of org sig: {2}".format(number_of_ensembles, T_per_ensemble, self.kernel.size))
                Ax = np.reshape(trim_signal, (int(size_of_stack),
                                int(self.kernel.size)), order='F')
                shape_Ax = Ax.shape
                shape_tHK = self.kernel.shape
                if shape_Ax[1] == shape_tHK[0]:
                    stack = np.matmul(Ax, self.kernel)  # create stack data
                    ensembles[:, index] = stack.T
                else:
                    print("fail stack, wrong size")
            elif index == (number_of_ensembles - 1):
                T_overlap = T_per_ensemble + overlap                   # get end overlap
                start_index = (index * (T_per_ensemble) - 2) * size_of_stack
                trim_signal = signal[start_index:]
                sub_signals.append(trim_signal)
                sub_samples.append(np.arange(start_index, signal.size))
                Ax = np.reshape(trim_signal, (int(size_of_stack),
                                int(self.kernel.size)), order='F')
                shape_Ax = Ax.shape
                shape_tHK = self.kernel.shape
                print("size of last: {0}, with overlap: {1}".format(start_index, signal.size))
                if shape_Ax[1] == shape_tHK[0]:
                    stack = np.matmul(Ax, self.kernel)  # create stack data
                    ensembles[:, index] = stack.T
                else:
                    print("fail stack, wrong size in last")
            else:
                T_overlap = T_per_ensemble + overlap
                start_index = (((index) * T_per_ensemble) - 2) * size_of_stack
                end_index = start_index + (T_overlap * size_of_stack)
                trim_signal = signal[start_index:end_index]
                sub_signals.append(trim_signal)
                sub_samples.append(np.arange(start_index, end_index))
                print("size of first: {0}, with overlap: {1}".format(start_index, signal.size))
                Ax = np.reshape(trim_signal, (int(size_of_stack),
                                int(self.kernel.size)), order='F')
                shape_Ax = Ax.shape
                shape_tHK = self.kernel.shape
                if shape_Ax[1] == shape_tHK[0]:
                    stack = np.matmul(Ax, self.kernel)  # create stack data
                    ensembles[:, index] = stack.T
                else:
                    print("fail stack, wrong size")
        ensembles[:, ::2] = ensembles[:, ::2] * -1   # make all stacks same polarity
        return ensembles, sub_signals, sub_samples


##################################################
# define methods
##################################################
# define methods
def getFFT(signal):
    """
       :rtype numpy array
       :return: fft of the filter kernal
    """
    v_fft = fftpack.fft(signal)
    return v_fft[0:int(v_fft.size / 2 - 1)] / np.max(v_fft)


def getFrequnceyResponse(signal):
    """
       :rtype numpy array
       :return: frequeny response of the filter kernal
    """
    v_fft = fftpack.fft(signal)
    amplitude = np.sqrt(v_fft.real**2 + v_fft.imag**2)
    return amplitude[(int(amplitude.size / 2 - 1)):] / np.max(amplitude)


def getCrossCorrelation(signal1, signal2):
    x_corr = np.correlate(signal1, signal2, mode='full')
    return x_corr / np.max(np.abs(x_corr))


def getCoherence(signal1, signal2, sample_rate):
    f, coherence = signal.coherence(signal1, signal2, sample_rate, nperseg=1024)
    return f, coherence


def getSpectralDensity(signal):
    """
       X_d(f)X_d^*(f)
    """
    signal_f = fftpack.fft(signal)
    signal_f_conj = np.conjugate(signal_f)
    S_xx = signal_f * signal_f_conj
    # return S_xx[0:int(S_xx.size / 2 - 1)] / np.max(np.abs(S_xx))
    return S_xx


def getPhaseResponse(signal):
    """
       :rtype numpy array
       :return: Rhase response of the filter kernal
    """
    v_fft = fftpack.fft(signal)
    phase = np.arctan2(v_fft.imag, v_fft.real)
    return phase[0:(phase.size / 2 - 1)]


def getCWT(sig, wavelet=None, widths=None):
    if wavelet is not None:
        cwtmatr = signal.cwt(sig, signal.ricker, widths)
    return cwtmatr


def getColeCole(mx_decay,
                init_cond,
                init_tau,
                step_factor,
                window_widths,
                window_weights):
    """
    takes in a decay and returns best fitting cole-cole
    note:
    A tool to calculate the time - domain voltage response[V(t) / Vo] for a
    homogenous Cole - Cole model of the earth using the digital linear filter
    formulation given by Guptasarma(Geophys, vol 47, pg 1575, 1982)
    """
    time = np.zeros(window_widths.size)  # initiates time
    # convert window widths to accumlative times specific to algorithm
    for i in range(time.size):
        time[i] = (np.sum(window_widths[0:i + 1]) / 2.0) / 1000.0
    c = np.zeros(12)                      # conductivity array
    v_cole = np.zeros(window_widths.size)  # best fit cole-cole
    tau = np.zeros(12)                    # time constant array
    err = np.zeros((12, 12))               # error matrix
    cole_m = np.zeros((12, 12))            # matrix of chargeabilities
    radius = step_factor                 # radius of array fill
    c[7] = init_cond                     # center of cond. array
    tau[7] = init_tau                    # center of tau array
    tau10 = np.log10(tau[7])             # log of time constant
    idx = np.arange(0, 12)
    c = c[7] + radius * (idx - 7) / 40.0  # fill cond. array
    tau = np.power(10.0,
                   (tau10 + radius * (idx - 7) /
                    2.))  # fill tau array
    # create filter
    areg = np.asarray([-3.82704, -3.56608, -3.30512, -3.04416,
                       -2.78320, -2.52224, -2.26128, -2.00032,
                       -1.73936, -1.47840, -1.21744, -.95648,
                       -.69552, -0.43456, -0.17360, 0.08736,
                       0.34832, 0.60928, 0.87024, 1.13120, 1.39216])
    # create 2nd filter
    preg = np.asarray([0.000349998, -0.000418371, 0.000772828,
                      -0.000171356, 0.001022172, 0.000897638,
                      0.002208974, 0.003844944, 0.006809040,
                      0.013029162, 0.022661391, 0.042972904,
                      0.075423603, 0.139346367, 0.234486236,
                      0.366178323, 0.284615486, -0.235691746,
                      0.046994188, -0.005901946, 0.000570165])
    v_cole = np.zeros(time.size)      # initiate decay array
    minErr = 0.5                      # signify initial Low error
    c_idx = 0                         # index of cond. of min err
    tau_idx = 0                       # index of cond. of min err

    # loop through the arrays of cond. and tau
    for i in range(c.size):
        for j in range(tau.size):
            ax = c[i] * np.pi / 2.0
            for win in range(mx_decay.size):
                v_temp = 0.0
                for n in range(areg.size):
                    w = np.power(10.0, (areg[n] - np.log10(time[win])))
                    ex = (w * tau[j])**c[i]
                    y = np.complex(ex * np.cos(ax), ex * np.sin(ax))
                    z = 1.0 - 1.0 / (1.0 + y)
                    # print(ex * np.cos(ax), np.real(z))
                    v_temp = v_temp + preg[n] * np.real(z)
                v_cole[win] = v_temp

            # calculate error
            norm_weights = np.sum(window_weights) / window_weights.size
            # serr = (np.sum(np.power((mx_decay - v_cole), 2) *
            #         window_weights) / norm_weights)
            err = np.sqrt(np.sum(np.power(np.abs(mx_decay - v_cole) *
                          window_weights, 2)) / window_weights.size)
            # np.sqrt(serr / (window_weights.size))
            if err < minErr:
                c_idx = i
                tau_idx = j
                minErr = (err)

            mx_conversion = 1000.                 # conversion for mV/V
            cole_m[i, j] = (np.sum(v_cole * window_widths) /
                            np.sum(window_widths)) * mx_conversion  # calcs Mx

            # go back and calculate best fit cole-cole curve and save it
            for win in range(mx_decay.size):
                v_temp = 0.0
                for n in range(areg.size):
                    w = np.power(10.0, (areg[n] - np.log10(time[win])))
                    ex = np.power(w * tau[tau_idx], c[c_idx])
                    y = np.complex(ex * np.cos(ax), ex * np.sin(ax))
                    z = 1.0 - 1.0 / (1.0 + y)
                    v_temp = v_temp + preg[n] * np.real(z)
                v_cole[win] = v_temp

    # calculate the percent diff
    # percent_diff = (
    #     np.mean((np.abs((v_cole - mx_decay) /
    #             ((v_cole + mx_decay) / 2.))) * window_weights))

    return c[c_idx], tau[tau_idx], cole_m[c_idx, tau_idx], minErr, v_cole


def padNextPower2(ts, num_ensembles=1):
    """
    Input: numpy time-series
    returns: a zeros padded time-series
    rtype: numpy array
    - if ensemble option is chosen, the time-series is arranged
      in a nxm numpy array

    NOTE: Primary meant for Qmag and MT processing
    """
    if num_ensembles == 1:
        next_2 = 1 << (ts.size - 1).bit_length()
        diff_in_size = next_2 - ts.size
        pad = np.zeros(diff_in_size)
        return np.concatenate((ts, pad), axis=0)
    elif num_ensembles > 1:
        return None
    else:
        return None


def getBestFitColeCole(mx_decay,
                       init_cond,
                       init_tau,
                       step_factor,
                       window_widths,
                       iterations,
                       timeCenter,
                       window_weights):
    c = init_cond
    tau = init_tau
    r = step_factor
    stored_error = 0
    min_iter = 3                               # min amount of iterations
    for iter in range(iterations):
        c, tau, M, error, vs = getColeCole(mx_decay,
                                           c,
                                           tau,
                                           r,
                                           timeCenter,
                                           window_weights)
        delta_error = abs(stored_error - error) / ((stored_error + error) / 2.)
        print("iter: %i | c: %f | tau: %f | M: %f | error: %f | delta: %f" %
              (iter, c, tau, M, error, delta_error))
        stored_error = error
        # r = r / 2.
        if delta_error > 0.002 or iter < min_iter:
            r = r / 2.
        elif delta_error < 0.002:
            print("convergence accomplished! DONE")
            break


def getWeightedVs(stack, window_start, window_end, attenuation):
    """
    takes in data and returns decay
    Input:
    stack = Half period stacked Voltage data
    window_start = an array of window start times
    window_end = an array of window end times

    """
    timebase = np.round(window_end[len(window_end) - 1] / 1000.0)
    timebase = timebase * 1000
    time = np.arange(0, timebase, (timebase / stack.size))
    vsDecay = np.zeros((len(window_end)))
    # loop through and do every window
    for win in range(len(window_start)):
        # find how many samples in first window
        cntr = 0
        # get time span for tail of windows
        start_tmp = window_start[win] - (
            0.45 * (window_end[win] - window_start[win]))
        if win == (vsDecay.size - 1):
            end_tmp = window_end[win]
        else:
            end_tmp = window_end[win] + (
                0.45 * (window_end[win] - window_start[win]))
        # print start_tmp, end_tmp
        for i in range(stack.size):
            if time[i] >= start_tmp and time[i] <= end_tmp:
                cntr += 1
        # create window wieghts
        indx1 = np.arange(0, attenuation)
        weights = 0.5 - (0.5 * np.cos((2 * np.pi * indx1) / (indx1.size - 1)))
        # create new weights
        Wg = np.zeros(cntr)
        start_Wg = (indx1.size / 2.0 - 1.0) - (cntr / 2.0) + 1
        for r in range(Wg.size):
            Wg[r] = weights[int(start_Wg) + r]
        # print Wg
        # create vector storing weighted values
        Vs_window = np.zeros(cntr)
        Vs_window_ave = np.zeros(cntr)
        # get window times
        w_idx = np.zeros(cntr)
        # assign total time and time step
        count = 0
        for i in range(time.size):
            if time[i] >= start_tmp and time[i] <= end_tmp:
                w_idx[count] = time[i]
                Vs_window[count] = stack[i] * -1 * Wg[count]
                Vs_window_ave[count] = stack[i] * -1
                count += 1
        sumWin = np.sum(Vs_window)
        # print Vs_window
        vsDecay[win] = sumWin / cntr

    return vsDecay


def createBruteStackWindow(num_points):
    num_points = num_points - 6
    tmp = np.ones(num_points)
    tmp = tmp * 4

    # create full filter kernal
    f1 = np.zeros(tmp.size + 4)
    f1[0] = 1
    f1[1] = 3
    f1[f1.size - 2] = 3
    f1[f1.size - 1] = 1

    for j in range(tmp.size):
        f1[j + 2] = tmp[j]

    f1 = f1 / (4.0 * (num_points - 2))

    return f1


def createMexHatWavelet(num_points, a):
    return signal.ricker(num_points, a)


def createSlepianWindow(num_size, attenuation):
    return signal.slepian(num_size, width=0.3)


def createKaiserWindow(num_taps, attenuation):
    """
    creates a Kaiser window
    Input:
    num_taps = number of taps for the requested window

    """
    return signal.kaiser(num_taps, beta=attenuation)


def mbessel(position, max_iter):
    seq_m = np.arange(0, max_iter)
    fact_m = factorial(seq_m, exact=False)
    summation = 0.0
    for i in range(max_iter):
        inc = np.power(1.0 / fact_m * np.power(position * 0.5, i), 2)
        frac = inc / summation
        summation += inc
        if frac < 0.001:
            break
    return summation


def cheby_poly(n, x):
    pos = 0.0
    if np.abs(x) <= 1.0:
        pos = np.cos(n * np.arccos(x))
    else:
        pos = np.cosh(n * np.arccosh(x))
    return pos


def createHanningWindow(num_points):
    """
    creates a Hanning window filter kernal
    Input:
    num_points = number of taps for the requested window

    """
    # num_points = num_points - 2
    indx1 = np.arange(1, num_points + 1).T          # create  sequence array
    filterWindow = 0.5 * (1 - np.cos((2 * np.pi /
                          (indx1.size - 1)) * indx1))  # creates window

    return filterWindow


def createChebyshevWindow(num_taps, attenuation):
    """
    creates a Chebyshev window
    Input:
    num_taps = number of taps for the requested window

    """
    # weights = np.zeros(num_taps)
    # summation = 0.0
    # max_value = 0.0
    # tg = np.power(10, (attenuation / 20))
    # x0 = np.cosh(1.0 / (num_taps - 1)) * np.arccosh(tg)
    # M = (num_taps - 1) / 2
    # if (num_taps % 2) == 0:
    #     M = M + 0.5
    # for nn in range(int(num_taps / 2 + 1)):
    #     n = nn - M
    #     summation = 0.0
    #     for i in range(int(M)):
    #         summation += cheby_poly(num_taps - 1.0,
    #                                 x0 * np.cos((2.0 *
    #                                             n * np.pi * i) / num_taps))
    #     weights[nn] = tg + 2 * summation
    #     weights[num_taps - nn - 1] = weights[nn]
    #     if weights[nn] > max_value:
    #         max_value = weights[nn]
    # for ii in range(num_taps):
    #     weights[ii] /= max_value
    return signal.chebwin(num_taps, at=attenuation)


def getPrimaryVoltage(start, end, stack):
        """
        Extracts the Vp of the signal.
        Input:
        start = percent of the on time to start calculation
        end = percent of on time to end calculation
        e.g 50% to 90%
        """
        sumStart = int((start / 100.0) * (stack.size / 2))  # start Vp calc
        sumEnd = int((end / 100.0) * (stack.size / 2))     # end of Vp calc
        Vp = np.sum(stack[sumStart:sumEnd]) / (sumEnd - sumStart)

        return Vp

def loadCsvDem(filepath):
    my_data = genfromtxt(filepath, delimiter=',', skip_header=1)
    utms = np.vstack((my_data[:, 2], my_data[:, 3]))
    elevations = my_data[:, 4] 
    grid = np.vstack((my_data[:, 0], my_data[:, 1]))
    return utms, grid, elevations

def plotCsvDemLocationDifferences(local=None, utm=None):
    if local is None or utm is None:
        print("[INFO] Please provide Locations of DEM!!!!")
    else:
        delta_gridx = local[0, 1:] - local[0, :-1]
        delta_gridy = local[1, 1:] - local[1, :-1]
        delta_utmx = utm[0, 1:] - utm[0, :-1]
        delta_utmy = utm[1, 1:] - utm[1, :-1]

        fig = plt.figure()
        ax1 = fig.add_subplot(221)
        ax2 = fig.add_subplot(222)
        ax3 = fig.add_subplot(223)
        ax4 = fig.add_subplot(224)

        ax1.plot(delta_gridx, '.')
        ax2.plot(delta_gridy, '.')
        ax3.plot(delta_utmx, '.')
        ax4.plot(delta_utmy, '.')

        ax1.set_title('Grid X difference')
        ax2.set_title('Grid Y difference')
        ax3.set_title('Easting difference')
        ax4.set_title('Northing difference')
        plt.show()

def plotLocalAndUtmGrids(local=None, utm=None):
    if local is None or utm is None:
        print("[INFO] Please provide locations info!!!")
    else:
        fig = plt.figure()
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)
        ax1.plot(local[0, :], local[1, :], '.')
        ax2.plot(utm[0, :], utm[1, :], '.')
        ax1.set_title('Local Grid Coords.')
        ax1.set_xlabel('X (m)')
        ax1.set_ylabel('Y (m)')
        ax2.set_xlabel('Easting (m)')
        ax2.set_ylabel('Northing (m)')
        ax2.set_title('UTM Coords.')
        plt.show()

def loadDias(fileName):
    """
    Function for loading a dias file and returns a
    "Patch" class complete with sources and recievers

    Input:
    fileName = complete path to data file

    """
    lines = 0
    text_file = open(fileName, "r")

    # determin how many lines in the file
    while text_file.readline():
            lines += 1
    text_file.close()

    # initiate a patch
    patch = Jpatch()
    # read header information
    text_file = open(fileName, "r")
    # initiate reading control variables
    currRdg = 0
    for i, line in enumerate(text_file):
        if i == 4:
            Varinfo = line.split()
            header5 = line
            # print(Varinfo)
        elif i == 0:
                header1 = line
        elif i == 1:
                header2 = line
                id_info = line.split()
                patch.assignPatchID(id_info[1])
        elif i == 2:
                header3 = line
        elif i == 3:
                header4 = line
        elif i > 4:
            # try:
                    datatxt = line.split()
                    # do some Jdatamanagment stuff
                    varFields = Jreadtxtline(Varinfo, datatxt)
                    # verify if line is a new reading
                    if varFields.RDG == currRdg:
                        # add the dipoles
                        # Idp = Jdata.JinDipole(varFields)
                        Vpdp = JvoltDipole(varFields)
                        Rdg.addVoltageDipole(Vpdp)
                    else:
                        # create a reading
                        Rdg = Jreading(varFields.RDG)
                        Idp = JinDipole(varFields)
                        Vpdp = JvoltDipole(varFields)
                        Rdg.addVoltageDipole(Vpdp)
                        Rdg.addInDipole(Idp)
                        # add reading to the patch
                        patch.addreading(Rdg)
                        currRdg = varFields.RDG
            # except:
                    # pass

    text_file.close()
    headers = [header1, header2, header3, header4, header5]
    patch.assignHeaderInfo(headers)
    return patch


# ===================================================
# Dias Data specific class
class JinDipole:
    """
    Class containing Source information

    Initiate with a structure containing location
    and Current value + error

    """

    def __init__(self, InDpInfo):
        self.Tx1File = InDpInfo.TxFile
        self.Tx1East = float(InDpInfo.Tx1East)
        self.Tx1North = float(InDpInfo.Tx1North)
        self.Tx1Elev = float(InDpInfo.Tx1Elev)
        self.Tx2East = float(InDpInfo.Tx2East)
        self.Tx2North = float(InDpInfo.Tx2North)
        self.Tx2Elev = float(InDpInfo.Tx2Elev)
        self.Tx1x = float(InDpInfo.Tx1x)
        self.Tx1y = float(InDpInfo.Tx1y)
        self.Tx2x = float(InDpInfo.Tx2x)
        self.Tx2y = float(InDpInfo.Tx2y)
        self.In = float(InDpInfo.In)
        self.In_err = float(InDpInfo.In_err)
        # print(InDpInfo.In_err)

    def getTxStack(self, stack_dir):
        rec_num = self.reading
        bin_file = stack_dir + "P1_R"+ rec_num + "_TX.raw"

        f = open(bin_file, "rb")          # open the file
        lines = f.read()                  # read data as string (can modify this to get header info
        return_char = 0.                  # beginning of binary data idex
        for idx in range(100):
            if lines[idx] == "\r":
                return_char = idx
        bin_start = return_char + 2       # increment 2 spaces to start of binary data
        data = []                         # initiate array for data
        with open(bin_file, 'rb') as f:   # open file for data extraction
            f.seek(bin_start,
                   os.SEEK_SET)          # seek to beginning of binary data
            while True:
                b = f.read(8)             # read 8 bytes at a time
                if not b:                 # break out if end of file
                    # eof
                    break
                data.append(
                    struct.unpack('d',
                                  b))  # store data

        return data


class JvoltDipole:
    """
    object containing voltage information

    """

    def __init__(self, VoltDpinfo):
        self.dipole = VoltDpinfo.DIPOLE
        self.reading = VoltDpinfo.RDG
        self.status = VoltDpinfo.Status
        self.Rx1x = float(VoltDpinfo.Rx1x)
        self.Rx1y = float(VoltDpinfo.Rx1y)
        self.Rx2x = float(VoltDpinfo.Rx2x)
        self.Rx2y = float(VoltDpinfo.Rx2y)
        self.Rx1File = VoltDpinfo.Rx1File
        self.Rx1East = float(VoltDpinfo.Rx1East)
        self.Rx1North = float(VoltDpinfo.Rx1North)
        self.Rx1Elev = float(VoltDpinfo.Rx1Elev)
        self.Rx2File = VoltDpinfo.Rx2File
        self.Rx2East = float(VoltDpinfo.Rx2East)
        self.Rx2North = float(VoltDpinfo.Rx2North)
        self.Rx2Elev = float(VoltDpinfo.Rx2Elev)
        self.coupling = float(VoltDpinfo.Coupling)
        self.K = float(VoltDpinfo.k)
        try:
            self.xcorr = float(VoltDpinfo.Xcorr)
        except:
            self.xcorr = -99.0
        try:
            self.eta = float(VoltDpinfo.M)
            self.c = float(VoltDpinfo.C)
            self.tau = float(VoltDpinfo.Tau)
        except:
            self.eta = -99.0
            self.c = -99.0
            self.tau = -99.0
        self.Sp = float(VoltDpinfo.Sp)
        try:
            self.Vp = float(VoltDpinfo.Vp)
        except:
            self.Vp = -999.9

        self.Vp_err = float(VoltDpinfo.Vp_err)
        self.Rho = float(VoltDpinfo.Rho)
        self.flagRho = VoltDpinfo.Rho_QC
        self.Stack = float(VoltDpinfo.Stack)
        try:
            self.Mx = float(VoltDpinfo.Mx)
        except:
            self.Mx = -99.9
        self.Mx_err = float(VoltDpinfo.Mx_err)
        self.flagMx = VoltDpinfo.Mx_QC
        self.flagBad = VoltDpinfo.Status
        self.TimeBase = VoltDpinfo.TimeBase
        self.Vs = np.asarray(VoltDpinfo.Vs)
        self.In = float(VoltDpinfo.In)
        self.In_err = float(VoltDpinfo.In_err)

    def getDipoleStack(self, stack_dir):
        node1_id = self.Rx1File[:2]
        node2_id = self.Rx2File[:2]
        rec_num = self.reading
        bin_file = stack_dir + "P1_R"+ rec_num + "_" + node1_id + "_" + node2_id + ".stk"

        f = open(bin_file, "rb")          # open the file
        lines = f.read()                  # read data as string (can modify this to get header info
        return_char = 0.                  # beginning of binary data idex
        for idx in range(100):
            if lines[idx] == "\r":
                return_char = idx
        bin_start = return_char + 2       # increment 2 spaces to start of binary data
        data = []                         # initiate array for data
        with open(bin_file, 'rb') as f:   # open file for data extraction
            f.seek(bin_start,
                   os.SEEK_SET)           # seek to beginning of binary data
            while True:
                b = f.read(8)             # read 8 bytes at a time
                if not b:                 # break out if end of file
                    # eof
                    break
                data.append(
                    struct.unpack('d',
                                  b))  # store data
                # if len(b) < 8:
                #     break
        return data

    def getXplotpoint(self, Idp, dipole_dipole=False):
        x = -9999.
        if dipole_dipole:
            mid_tx = (Idp.tx1x + Idp.tx2x) / 2.0
            if (self.Rx1x > mid_tx):
                x = mid_tx + ((self.Rx1x - mid_tx) / 2.0)
            elif (self.Rx1x < mid_tx):
                x = Idp.Tx1x - ((Idp.Tx1x - self.Rx1x) / 2.0)
        else:
            if (self.Rx1x > Idp.Tx1x):
                x = Idp.Tx1x + ((self.Rx1x - Idp.Tx1x) / 2.0)
            elif (self.Rx1x < Idp.Tx1x):
                x = Idp.Tx1x - ((Idp.Tx1x - self.Rx1x) / 2.0)
        return x

    def getYplotpoint(self, Idp, dipole_dipole=False):
        y = -9999.
        if dipole_dipole:
            mid_tx = (Idp.tx1y + Idp.tx2y) / 2.0
            if (self.Rx1y > Idp.Tx1y):
                y = mid_tx + ((self.Rx1y - mid_tx) / 2.0)
            elif (self.Rx1y < mid_tx):
                y = mid_tx - ((mid_tx - self.Rx1y) / 2.0)
        else:
            if (self.Rx1y > Idp.Tx1y):
                y = Idp.Tx1y + ((self.Rx1y - Idp.Tx1y) / 2.0)
            elif (self.Rx1y < Idp.Tx1y):
                y = Idp.Tx1y - ((Idp.Tx1y - self.Rx1y) / 2.0)
        return y

    def getZplotpoint(self, Idp, dipole_dipole=False):
        z = -9999.
        if dipole_dipole:
            mid_tx = (Idp.tx1x + Idp.tx2x) / 2.0
            z = -(abs(mid_tx - self.Rx1x)) / 2.0
        else:
            r = np.sqrt((Idp.Tx1x - self.Rx2x)**2 +
                        (Idp.Tx1y - self.Rx2y)**2 +
                        (Idp.Tx1Elev - self.Rx2Elev)**2)
            z = -(r / 3.)
        return z

    def calcSensitivity(self, Idp, ne, sw, dx,
                        sw_search,
                        ne_search):
        xa = Idp.Tx2East
        ya = Idp.Tx2North
        za = Idp.Tx2Elev
        xb = Idp.Tx1East
        yb = Idp.Tx1North
        zb = Idp.Tx1Elev
        xm = self.Rx1East
        ym = self.Rx1North
        zm = self.Rx1Elev
        xn = self.Rx2East
        yn = self.Rx2North
        zn = self.Rx2Elev
        X = np.arange(sw[0], ne[0], dx)
        Y = np.arange(sw[1], ne[1], dx)
        Z = np.arange(sw[2], ne[2], dx)
        x, y, z = np.meshgrid(X, Y, Z)
        # compactify the calculation
        arg1 = (((x - xa) * (x - xm) + (y - ya) * (y - ym) +
                (z - za) * (z - zm)) *
                (((x - xa)**2 + (y - ya)**2 + (z - za)**2)**-1.5) *
                (((x - xm)**2 + (y - ym)**2 + (z - zm)**2)**-1.5))
        arg2 = (((x - xa) * (x - xn) + (y - ya) * (y - yn) +
                (z - za) * (z - zn)) *
                (((x - xa)**2 + (y - ya)**2 + (z - za)**2)**-1.5) *
                (((x - xn)**2 + (y - yn)**2 + (z - zn)**2)**-1.5))
        arg3 = (((x - xb) * (x - xm) + (y - yb) * (y - ym) +
                (z - zb) * (z - zm)) *
                (((x - xb)**2 + (y - yb)**2 + (z - zb)**2)**-1.5) *
                (((x - xm)**2 + (y - ym)**2 + (z - zm)**2)**-1.5))
        arg4 = (((x - xb) * (x - xn) + (y - yb) * (y - yn) +
                (z - zb) * (z - zn)) *
                (((x - xb)**2 + (y - yb)**2 + (z - zb)**2)**-1.5) *
                (((x - xn)**2 + (y - yn)**2 + (z - zn)**2)**-1.5))
        # evaluate sensitivity equation
        J = 1 / (4 * np.pi**2) * (arg1 - arg2 - arg3 + arg4)
        # get the span of the most sensitive layer
        max_J = np.max(np.abs(J)) * 0.91
        levels = J.shape[2]
        slices_x = np.zeros((levels, 2))
        slices_y = np.zeros((levels, 2))
        for l in range(levels):
            z = J[:, :, l]
            a = np.abs(z) > (np.max(np.abs(z)) * 0.61)
            X_ = x[a]
            Y_ = y[a]
            try:
                slices_x[l, :] = np.asarray([np.max(X_), np.min(X_)])
                slices_y[l, :] = np.asarray([np.max(Y_), np.min(Y_)])
            except ValueError:
                slices_x[l, :] = slices_x[0, :]
                slices_y[l, :] = slices_y[0, :]

        sbx = [np.max(slices_x[:, 0]), np.min(slices_x[:, 1])]
        sby = [np.max(slices_y[:, 0]), np.min(slices_y[:, 1])]
        # determine if it overlaps with the area in question
        true_x = False
        true_y = False
        x_s = sw_search[0]
        x1_s = ne_search[0]
        y_s = sw_search[1]
        y1_s = ne_search[1]
        if (sbx[1] < x_s < sbx[0]):
            true_x = True
        if (sbx[1] < x1_s < sbx[0]):
            true_x = True
        if (sby[1] < y_s < sby[0]):
            true_y = True
        if (sby[1] < y1_s < sby[0]):
            true_y = True
        # print("Max sensitive x comp. {0} - min: {1}".format(sbx[0], sbx[1]))
        # print("Max search x comp. {0} - min: {1}".format(x1_s, x_s))
        # print("Max sensitive y comp. {0} - min: {1}".format(sby[0], sby[1]))
        # print("Max search y comp. {0} - min: {1}".format(y1_s, y_s))
        # print(true_x, true_y)
        if true_x and true_y:
            print("found sensitive dipole")
            self.flagRho = "Reject"
            self.flagMx = "Reject"

        return J

    def getDepthPoint2D(self, Idp, direction='x'):
        z = np.nan
        if direction is 'x':
            z = -(abs(Idp.Tx1x - self.Rx1x)) / 2.0
        elif direction is 'y':
            z = -(abs(Idp.Tx1y - self.Rx1y)) / 2.0
        return z

    def getNlevel(self, Idp, a, direction='x'):
        n = np.nan
        if direction is 'x':
            n = (abs(Idp.Tx1x - self.Rx1x)) / a
        elif direction is 'y':
            n = -(abs(Idp.Tx1y - self.Rx1y)) / a
        return n

    def getAseperation(self):
        r1 = ((self.Rx1East - self.Rx2East)**2 +
              (self.Rx1North - self.Rx2North)**2 +
              (self.Rx1Elev - self.Rx2Elev)**2)**0.5
        return r1

    def getAseperationLocal(self, direction='x'):
        a = np.nan
        if direction is 'x':
            a = np.abs(self.Rx2x - self.Rx1x)
        elif direction is 'y':
            a = np.abs(self.Rx2y - self.Rx1y)
        return a

    def calcGeoFactor(self, Idp):
        r1 = ((self.Rx1East - Idp.Tx1East)**2 +
              (self.Rx1North - Idp.Tx1North)**2 +
              (self.Rx1Elev - Idp.Tx1Elev)**2)**0.5
        r2 = ((self.Rx2East - Idp.Tx1East)**2 +
              (self.Rx2North - Idp.Tx1North)**2 +
              (self.Rx2Elev - Idp.Tx1Elev)**2)**0.5
        r3 = ((self.Rx1East - Idp.Tx2East)**2 +
              (self.Rx1North - Idp.Tx2North)**2 +
              (self.Rx1Elev - Idp.Tx2Elev)**2)**0.5
        r4 = ((self.Rx2East - Idp.Tx2East)**2 +
              (self.Rx2North - Idp.Tx2North)**2 +
              (self.Rx2Elev - Idp.Tx2Elev)**2)**0.5
        try:
            gf = 1 / ((1 / r1 - 1 / r2) - (1 / r3 - 1 / r4))
        except ZeroDivisionError:
            gf = 0.0
        return 2 * np.pi * gf

    def calcRho(self, Idp):
        r1 = ((self.Rx1East - Idp.Tx1East)**2 +
              (self.Rx1North - Idp.Tx1North)**2 +
              (self.Rx1Elev - Idp.Tx1Elev)**2)**0.5
        r2 = ((self.Rx2East - Idp.Tx1East)**2 +
              (self.Rx2North - Idp.Tx1North)**2 +
              (self.Rx2Elev - Idp.Tx1Elev)**2)**0.5
        r3 = ((self.Rx1East - Idp.Tx2East)**2 +
              (self.Rx1North - Idp.Tx2North)**2 +
              (self.Rx1Elev - Idp.Tx2Elev)**2)**0.5
        r4 = ((self.Rx2East - Idp.Tx2East)**2 +
              (self.Rx2North - Idp.Tx2North)**2 +
              (self.Rx2Elev - Idp.Tx2Elev)**2)**0.5
        try:
            gf = 1 / ((1 / r1 - 1 / r2) - (1 / r3 - 1 / r4))
        except ZeroDivisionError:
            gf = 0.0
        Vp = self.Vp
        # print("Vp: {0}".format(self.Vp))
        rho = (Vp / self.In) * 2 * np.pi * gf
        # self.Rho = rho
        return rho

    def calcMx(self, window_widths, start, stop):
        widths = window_widths[start:stop]
        vs_ = self.Vs[start:stop]
        mx_calc = np.sum(vs_ * widths) / np.sum(widths)
        return mx_calc / (self.Vp / 1000.0) 


class Jreading:
    """
    Class to handle current and voltage dipole
    information for a given source

    """
    def __init__(self, mem):
        self.MemNumber = mem
        self.Vdp = []
        self.win_width = []
        self.win_start = []
        self.win_end = []
        self.node_db = []
        self.node_locs = []
        self.node_ids = []
        self.angles = []

    # method for adding voltage dipoles
    def addVoltageDipole(self, JVdp):
        self.Vdp.append(JVdp)

    # method for assigning Current dipole
    def addInDipole(self, JIdp):
        self.Idp = JIdp

    # method for creating a node database
    def createNodeDB(self):
        nodes = []
        locs = []
        node_num = []
        for dp in range(len(self.Vdp)):
            node1_id = self.Vdp[dp].Rx1File[:2]
            node2_id = self.Vdp[dp].Rx2File[:2]
            if not node1_id in nodes:
                nodes.append(node1_id)
                pos = np.asarray([self.Vdp[dp].Rx1East, self.Vdp[dp].Rx1North])
                locs.append(pos)
                node_num.append(self.Vdp[dp].Rx1File)
            if not node2_id in nodes:
                nodes.append(node2_id)
                pos = np.asarray([self.Vdp[dp].Rx2East, self.Vdp[dp].Rx2North])
                locs.append(pos)
                node_num.append(self.Vdp[dp].Rx2File)
        # assign node database
        self.node_db = nodes
        self.node_locs = locs
        self.node_ids = node_num

    # method for getting angles of possible dipole given a spacing
    def calcAngles(self, min_dipole, max_dipole):
        # dipoles = [] len(self.node_locs
        if len(self.node_locs) > 1:
            points = np.asarray(self.node_locs)
            for idx in range(len(self.node_locs)):
                # calculate difference to all points
                diff = self.node_locs[idx] - points
                # calculate distances
                dists = (np.sum(diff**2, axis=1))**0.5
                for index in range(dists.size):
                    if min_dipole < dists[index] < max_dipole:
                        # get components differences
                        x_diff = diff[index, 0]
                        y_diff = diff[index, 1]
                        # calc angle from north
                        theta = (np.arctan2(x_diff, y_diff)) * 180 / np.pi
                        if theta < 0:
                            theta = 360 + theta
                        self.angles.append(theta)
        else:
            print("node locations have not been initiated")


class Jpatch:
    """
    Class to hold source information for a given data patch

    """

    def __init__(self):
        self.readings = []

    def addreading(self, Jrdg):
        self.readings.append(Jrdg)

    def assignPatchID(self, id):
        self.ID = id

    def assignHeaderInfo(self, headerLines):
        """
        Method for processing the header lines
        of a Dias data file. (e.g IP times, project name, etc.)

        Input: an array of header lines from file

        """

        self.headers = headerLines         # assigns headers to patch class
        # process IP times from header
        timing_string = self.headers[2].split(' ')[1]
        timings = timing_string.split(';')[0]
        timing = timings.split(',')
        self.window_start = np.zeros(len(timing))
        self.window_end = np.zeros(len(timing))
        self.window_center = np.zeros(len(timing))
        self.window_width = np.zeros(len(timing))
        # loops each available IP window and calculate width and centre
        for i in range(len(timing)):
            wintime = timing[i].split(':')
            self.window_start[i] = float(wintime[0])
            self.window_end[i] = float(wintime[1])
            self.window_center[i] = (self.window_start[i] +
                                     self.window_end[i]) / 2.0
            self.window_width[i] = (self.window_end[i] -
                                    self.window_start[i])

    def assignWindowWidthsToDipoles(self):
        """
        assigns widths of windows of Vs decay to each dipole

        """
        num_rdg = len(self.readings)
        for k in range(num_rdg):
            self.readings[k].win_width = self.window_width
            self.readings[k].win_start = self.window_start
            self.readings[k].win_end = self.window_end

    def rejectByVsError(self):
        num_rdg = len(self.readings)
        for rdg in range(num_rdg):
            num_dipole = len(self.readings[rdg].Vdp)
            for dp in range(num_dipole):
                if self.readings[rdg].Vdp[dp].Mx_err == -99.9:
                        self.readings[rdg].Vdp[dp].flagMx = "Reject"

    def get3Dangles(self, bin_width=20, dipole_max=None, dipole_min=None):
        dipole_angles = []
        theta = []
        cnt = 0
        if dipole_max is not None:
            num_rdg = len(self.readings)
            for rdg in range(num_rdg):
                # get a node database for a reading
                self.readings[rdg].createNodeDB()
                self.readings[rdg].calcAngles(dipole_min, dipole_max)
                dipole_angles.append(self.readings[rdg].angles)
                for idx in range(len(dipole_angles[cnt])):
                    theta.append(dipole_angles[cnt][idx])
                cnt = cnt + 1
        else:
            print("[OUTPUT] - dipole sizes not set")
        return theta

    # method for making node database for entire project
    def getNodeDataBase(self):
        num_rdg = len(self.readings)
        nodes = []
        # locs = []
        for rdg in range(num_rdg):
            self.readings[rdg].createNodeDB()
            num_dipole = len(self.readings[rdg].Vdp)
            for dp in range(num_dipole):
                node1_id = self.readings[rdg].Vdp[dp].Rx1File[:2]
                node2_id = self.readings[rdg].Vdp[dp].Rx2File[:2]
                if not node1_id in nodes:
                    nodes.append(node1_id)
                    pos = np.asarray([self.readings[rdg].Vdp[dp].Rx1East,
                                     self.readings[rdg].Vdp[dp].Rx1North])
                    # locs.append(pos)
                if not node2_id in nodes:
                    nodes.append(node2_id)
                    pos = np.asarray([self.readings[rdg].Vdp[dp].Rx2East,
                                     self.readings[rdg].Vdp[dp].Rx2North])
                    # locs.append(pos)
        return nodes

    def checkForMissingRecords(self, record_log=None):
        if record_log is not None:
            rec_out_file = record_log + ".rlog"
            text_file = open(record_log, "r")
            outfile = open(rec_out_file, "w+")
            # retrieve and assign the data
            rec_log = []
            current_rec = -1
            for line in text_file:
                spl = line.split(':')
                if int(spl[1]) != current_rec:
                    rec_log.append(int(spl[1]))
                    current_rec = int(spl[1])
            for rdg in range(len(self.readings)):
                found_match = False
                for rec in rec_log:
                    if int(self.readings[rdg].MemNumber) == rec:
                        found_match = True
                if not found_match:
                    outfile.write('rec %s found no match\n' % self.readings[rdg].MemNumber)
            print("[OUTPUT]  Log file written to same file as records log location")
        else:
            print(      "Please supply a records log file!!!!!")

    def checkContinuityInNodeFiles(self, path=None):
        if path is not None:
            nodes = self.getNodeDataBase()
            nodes.sort()
            # print(nodes)
            outfile = open(path, "w+")
            for idx in range(len(nodes)):
                nid = nodes[idx]
                files = []
                hold_num = -99
                hold_id = ""
                hold_mem = -99
                for rdg in range(len(self.readings)):
                    for ids in range(len(self.readings[rdg].node_ids)):
                        if nid == self.readings[rdg].node_ids[ids][:2]:
                            if hold_num == -99:
                                outfile.write('node %s mem start: %s\n' % (nid, self.readings[rdg].MemNumber))
                                hold_num = int(self.readings[rdg].node_ids[ids][2:8])
                                hold_mem = int(self.readings[rdg].MemNumber)
                                hold_id = self.readings[rdg].node_ids[ids]
                            else:
                                diff = int(self.readings[rdg].node_ids[ids][2:8]) - hold_num
                                if diff > 1:
                                    outfile.write('found gap: %s mem: %i & %s mem: %s\n' % (hold_id,
                                                   hold_mem, self.readings[rdg].node_ids[ids],
                                                   self.readings[rdg].MemNumber))
                                    hold_num = int(self.readings[rdg].node_ids[ids][2:8])
                                    hold_mem = int(self.readings[rdg].MemNumber)
                                    hold_id = self.readings[rdg].node_ids[ids]
            # write node database collected from db file
            for ids in nodes:
                outfile.write('%s\n' % ids)
            print("It be done!, file written")
        else:
            print("please supply path to database file")

    def rejectDataFromSensitiveArea(self, dx, sw_s, ne_s, outFile=None):
        """
           This function takes in an area of concern and removes data
           sensitive to the defined area
        """
        if outFile is not None:
            rx_locs = self.getDipoles()
            tx_locs = self.getSources()
            sw = np.array([np.min(rx_locs[:, 0]), np.min(tx_locs[:, 1]),
                  np.min(rx_locs[:, 2]) - 200])
            ne = np.array([np.max(rx_locs[:, 0]), np.max(rx_locs[:, 1]),
                  np.max(rx_locs[:, 2])])
            X = np.arange(sw[0], ne[0], dx)
            Y = np.arange(sw[1], ne[1], dx)
            Z = np.arange(sw[2], ne[2], dx)
            J_ = np.zeros((Y.size, X.size, Z.size))
            for rdg in range(len(self.readings)):
                current = self.readings[rdg].Idp
                print("rec num: {0}".format(self.readings[rdg].MemNumber))
                for dp in range(len(self.readings[rdg].Vdp)):
                    J_ += self.readings[rdg].Vdp[dp].calcSensitivity(current,ne, sw, dx, sw_search=sw_s, ne_search=ne_s)
            
            start_time = self.window_start[0]
            end_time = self.window_end[0]
            start_inds = (window_start >= start_time)
            stop_inds = (window_end <= end_time)
            strt_time = window_start[start_inds]
            stop_time = window_end[stop_inds]
            # write the data to file
            self.writeColeColeSEDat(outname, strt_time[0],
                                     stop_time[stop_time.size - 1])
            print("New File writte with sensitive area removed")
        else:
            print("     Please supply an output file path!")


    def writeColeColeSEDat(self, outname, start_time, end_time):
        """
        Writes Cole-Cole data from stretched exponential to file for DiasQC
        """
        out_file = open(outname, "w+")
        format_jup = '%6s %8s %8s %8i %13s %12.3f %12.3f %12.3f %12.0f %12.0f %12.3f %12.3f %12.3f %12.0f %12.0f %14s %12.3f %12.3f %12.3f %12.0f %12.0f %14s %12.3f %12.3f %12.3f %12.0f %12.0f %14.1f %10.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12s %14.3f %14.3f %11.0f %13s %8.3f %8.3f %8s %12.3e %12.3e %10.2f '
        # write the headers
        for i in range(len(self.headers)):
            if i == 2:
                broken = self.headers[2].split()
                self.headers[i] = broken[0] + " " + broken[1] + " - " + str(int(start_time)) + ":" + str(int(end_time)) + "(ms)\n"
                out_file.write('%s' % self.headers[i])
            elif i < (len(self.headers) - 1):
                out_file.write('%s' % self.headers[i])
            else:
                # fix headers
                hdr_temp = self.headers[i].split("Vs01")
                self.headers[i] = hdr_temp[0] + "   C          Tau          M         Vs01" + hdr_temp[1]
                out_file.write('%s' % self.headers[i])

        for rec in range(len(self.readings)):
            num_dipole = len(self.readings[rec].Vdp)
            for dp in range(num_dipole):
                if self.readings[rec].Vdp[dp].c is None:
                    cole_con = -99.9
                    cole_tau = -99.9
                    cole_eta = -99.9
                else:
                    cole_con = self.readings[rec].Vdp[dp].c
                    cole_tau = self.readings[rec].Vdp[dp].tau
                    cole_eta = self.readings[rec].Vdp[dp].eta
                if np.abs(self.readings[rec].Vdp[dp].Mx_err) > 999.0:
                    self.readings[rec].Vdp[dp].Mx_err = -99.9
                    self.readings[rec].Vdp[dp].flagMx = "Reject"
                if np.abs(self.readings[rec].Vdp[dp].Mx) > 999.0:
                    self.readings[rec].Vdp[dp].Mx = -99.9
                    self.readings[rec].Vdp[dp].flagMx = "Reject"
                out_file.write(format_jup % (str(self.readings[rec].Vdp[dp].reading),
                                             str(self.readings[rec].Vdp[dp].dipole),
                                             str(self.readings[rec].Vdp[dp].status),
                                             int(self.readings[rec].Vdp[dp].getAseperation()),
                                             self.readings[rec].Idp.Tx1File,
                                             self.readings[rec].Idp.Tx1East,
                                             self.readings[rec].Idp.Tx1North,
                                             self.readings[rec].Idp.Tx1Elev,
                                             self.readings[rec].Idp.Tx1x,
                                             self.readings[rec].Idp.Tx1y,
                                             self.readings[rec].Idp.Tx2East,
                                             self.readings[rec].Idp.Tx2North,
                                             self.readings[rec].Idp.Tx2Elev,
                                             self.readings[rec].Idp.Tx2x,
                                             self.readings[rec].Idp.Tx2y,
                                             self.readings[rec].Vdp[dp].Rx1File,
                                             self.readings[rec].Vdp[dp].Rx1East,
                                             self.readings[rec].Vdp[dp].Rx1North,
                                             self.readings[rec].Vdp[dp].Rx1Elev,
                                             self.readings[rec].Vdp[dp].Rx1x,
                                             self.readings[rec].Vdp[dp].Rx1y,
                                             self.readings[rec].Vdp[dp].Rx2File,
                                             self.readings[rec].Vdp[dp].Rx2East,
                                             self.readings[rec].Vdp[dp].Rx2North,
                                             self.readings[rec].Vdp[dp].Rx2Elev,
                                             self.readings[rec].Vdp[dp].Rx2x,
                                             self.readings[rec].Vdp[dp].Rx2y,
                                             0.0,
                                             self.readings[rec].Vdp[dp].xcorr,
                                             self.readings[rec].Vdp[dp].Sp,
                                             self.readings[rec].Vdp[dp].Vp,
                                             self.readings[rec].Vdp[dp].Vp_err,
                                             self.readings[rec].Vdp[dp].In,
                                             self.readings[rec].Vdp[dp].In_err,
                                             self.readings[rec].Vdp[dp].calcRho(self.readings[rec].Idp),
                                             self.readings[rec].Vdp[dp].flagRho,
                                             self.readings[rec].Vdp[dp].calcGeoFactor(self.readings[rec].Idp),
                                             self.readings[rec].Vdp[dp].coupling,
                                             self.readings[rec].Vdp[dp].Stack,
                                             self.readings[rec].Vdp[dp].TimeBase,
                                             self.readings[rec].Vdp[dp].Mx,
                                             self.readings[rec].Vdp[dp].Mx_err,
                                             self.readings[rec].Vdp[dp].flagMx,
                                             cole_con,
                                             cole_tau,
                                             cole_eta))
                for win in range(self.readings[rec].Vdp[dp].Vs.size):
                    if win < (self.readings[rec].Vdp[dp].Vs.size - 1):
                        out_file.write('%12.4f' % (self.readings[rec].Vdp[dp].Vs[win]))
                    else:
                        out_file.write('%12.4f\n' % (self.readings[rec].Vdp[dp].Vs[win]))

    def writeColeColeDat(self, eta, tau, c,
                         error, indicies, outname, start_time, end_time):
        """
        Writes Cole-Cole data from stretched exponential to file for DiasQC
        """
        for idx in range(indicies.shape[0]):
            rdg = indicies[idx, 0]
            dp = indicies[idx, 1]
            self.readings[rdg].Vdp[dp].c = c[idx]
            self.readings[rdg].Vdp[dp].eta = eta[idx]
            self.readings[rdg].Vdp[dp].tau = tau[idx]
            if error[idx] == -99.9:
                self.readings[rdg].Vdp[dp].Mx_err = error[idx]
                self.readings[rdg].Vdp[dp].flagMx = "Reject"
            else:
                self.readings[rdg].Vdp[dp].Mx_err = error[idx]
                self.readings[rdg].Vdp[dp].flagMx = "Accept"

        out_file = open(outname, "w+")
        format_jup = '%6s %8s %8s %8i %13s %12.3f %12.3f %12.3f %12.0f %12.0f %12.3f %12.3f %12.3f %12.0f %12.0f %14s %12.3f %12.3f %12.3f %12.0f %12.0f %14s %12.3f %12.3f %12.3f %12.0f %12.0f %14.1f %10.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12s %14.3f %14.3f %11.0f %13s %8.3f %8.3f %8s %12.3e %12.3e %10.2f '
        # write the headers
        for i in range(len(self.headers)):
            if i == 2:
                broken = self.headers[2].split()
                self.headers[i] = broken[0] + " " + broken[1] + " - " + str(int(start_time)) + ":" + str(int(end_time)) + "(ms)\n"
                out_file.write('%s' % self.headers[i])
            elif i < (len(self.headers) - 1):
                out_file.write('%s' % self.headers[i])
            else:
                # fix headers
                hdr_temp = self.headers[i].split("Vs01")
                self.headers[i] = hdr_temp[0] + "   C          Tau          M         Vs01" + hdr_temp[1]
                out_file.write('%s' % self.headers[i])

        for rec in range(len(self.readings)):
            num_dipole = len(self.readings[rec].Vdp)
            for dp in range(num_dipole):
                if self.readings[rec].Vdp[dp].c is None:
                    cole_con = -99.9
                    cole_tau = -99.9
                    cole_eta = -99.9
                else:
                    cole_con = self.readings[rec].Vdp[dp].c
                    cole_tau = self.readings[rec].Vdp[dp].tau
                    cole_eta = self.readings[rec].Vdp[dp].eta
                if np.abs(self.readings[rec].Vdp[dp].Mx_err) > 999.0:
                    self.readings[rec].Vdp[dp].Mx_err = -99.9
                    self.readings[rec].Vdp[dp].flagMx = "Reject"
                if np.abs(self.readings[rec].Vdp[dp].Mx) > 999.0:
                    self.readings[rec].Vdp[dp].Mx = -99.9
                    self.readings[rec].Vdp[dp].flagMx = "Reject"
                out_file.write(format_jup % (str(self.readings[rec].Vdp[dp].reading),
                                             str(self.readings[rec].Vdp[dp].dipole),
                                             str(self.readings[rec].Vdp[dp].status),
                                             int(self.readings[rec].Vdp[dp].getAseperation()),
                                             self.readings[rec].Idp.Tx1File,
                                             self.readings[rec].Idp.Tx1East,
                                             self.readings[rec].Idp.Tx1North,
                                             self.readings[rec].Idp.Tx1Elev,
                                             self.readings[rec].Idp.Tx1x,
                                             self.readings[rec].Idp.Tx1y,
                                             self.readings[rec].Idp.Tx2East,
                                             self.readings[rec].Idp.Tx2North,
                                             self.readings[rec].Idp.Tx2Elev,
                                             self.readings[rec].Idp.Tx2x,
                                             self.readings[rec].Idp.Tx2y,
                                             self.readings[rec].Vdp[dp].Rx1File,
                                             self.readings[rec].Vdp[dp].Rx1East,
                                             self.readings[rec].Vdp[dp].Rx1North,
                                             self.readings[rec].Vdp[dp].Rx1Elev,
                                             self.readings[rec].Vdp[dp].Rx1x,
                                             self.readings[rec].Vdp[dp].Rx1y,
                                             self.readings[rec].Vdp[dp].Rx2File,
                                             self.readings[rec].Vdp[dp].Rx2East,
                                             self.readings[rec].Vdp[dp].Rx2North,
                                             self.readings[rec].Vdp[dp].Rx2Elev,
                                             self.readings[rec].Vdp[dp].Rx2x,
                                             self.readings[rec].Vdp[dp].Rx2y,
                                             0.0,
                                             self.readings[rec].Vdp[dp].xcorr,
                                             self.readings[rec].Vdp[dp].Sp,
                                             self.readings[rec].Vdp[dp].Vp,
                                             self.readings[rec].Vdp[dp].Vp_err,
                                             self.readings[rec].Vdp[dp].In,
                                             self.readings[rec].Vdp[dp].In_err,
                                             self.readings[rec].Vdp[dp].calcRho(self.readings[rec].Idp),
                                             self.readings[rec].Vdp[dp].flagRho,
                                             self.readings[rec].Vdp[dp].calcGeoFactor(self.readings[rec].Idp),
                                             self.readings[rec].Vdp[dp].coupling,
                                             self.readings[rec].Vdp[dp].Stack,
                                             self.readings[rec].Vdp[dp].TimeBase,
                                             self.readings[rec].Vdp[dp].Mx,
                                             self.readings[rec].Vdp[dp].Mx_err,
                                             self.readings[rec].Vdp[dp].flagMx,
                                             cole_con,
                                             cole_tau,
                                             cole_eta))
                for win in range(self.readings[rec].Vdp[dp].Vs.size):
                    if win < (self.readings[rec].Vdp[dp].Vs.size - 1):
                        out_file.write('%12.4f' % (self.readings[rec].Vdp[dp].Vs[win]))
                    else:
                        out_file.write('%12.4f\n' % (self.readings[rec].Vdp[dp].Vs[win]))

    def getApparentResistivity(self, reject="app_rho", calc=True):
        """
        Exports all the apparent resistivity data

        Output:
        numpy array [value]

        """
        resistivity_list = []
        num_rdg = len(self.readings)
        for k in range(num_rdg):
            num_dipole = len(self.readings[k].Vdp)
            for j in range(num_dipole):
                if reject is "app_rho":
                    if self.readings[k].Vdp[j].flagRho == "Accept":
                        if calc:
                            rho_a = self.readings[k].Vdp[j].calcRho(self.readings[k].Idp)
                            resistivity_list.append(rho_a)
                        else:
                            rho_a = self.readings[k].Vdp[j].Rho
                            resistivity_list.append(rho_a)
                elif reject is "app_mx":
                    if self.readings[k].Vdp[j].flagMx == "Accept":
                        if calc:
                            rho_a = self.readings[k].Vdp[j].calcRho(self.readings[k].Idp)
                            resistivity_list.append(rho_a)
                        else:
                            rho_a = self.readings[k].Vdp[j].Rho
                            resistivity_list.append(rho_a)

        return np.vstack(resistivity_list)

    def getFrequencyComponent(self, reject="app_rho"):
        """
        Exports all the frequency component data

        Output:
        numpy array [value]

        """
        c_list = []
        num_rdg = len(self.readings)
        for k in range(num_rdg):
            num_dipole = len(self.readings[k].Vdp)
            for j in range(num_dipole):
                if reject is "app_rho":
                    if self.readings[k].Vdp[j].flagRho == "Accept":
                        c_list.append(self.readings[k].Vdp[j].c)
                elif reject is "app_mx":
                    if self.readings[k].Vdp[j].flagMx == "Accept":
                        c_list.append(self.readings[k].Vdp[j].c)

        return np.vstack(c_list)

    def getTimeConstants(self, reject="app_rho"):
        """
        Exports all the time constant data

        Output:
        numpy array [value]

        """
        tau_list = []
        num_rdg = len(self.readings)
        for k in range(num_rdg):
            num_dipole = len(self.readings[k].Vdp)
            for j in range(num_dipole):
                if reject is "app_rho":
                    if self.readings[k].Vdp[j].flagRho == "Accept":
                        tau_list.append(self.readings[k].Vdp[j].tau)
                elif reject is "app_mx":
                    if self.readings[k].Vdp[j].flagMx == "Accept":
                        tau_list.append(self.readings[k].Vdp[j].tau)

        return np.vstack(tau_list)

    def getColeColeMx(self, reject="Rho"):
        """
        Exports all the eta data

        Output:
        numpy array [value]

        """
        eta_list = []
        num_rdg = len(self.readings)
        for k in range(num_rdg):
            num_dipole = len(self.readings[k].Vdp)
            for j in range(num_dipole):
                if reject is "app_rho":
                    if self.readings[k].Vdp[j].flagRho == "Accept":
                        eta_list.append(self.readings[k].Vdp[j].eta)
                elif reject is "app_mx":
                    if self.readings[k].Vdp[j].flagMx == "Accept":
                        eta_list.append(self.readings[k].Vdp[j].eta)

        return np.vstack(eta_list)

    def getVsErrors(self, reject=None):
        """
        Exports all the eta data

        Output:
        numpy array [value]

        """
        error_list = []
        num_rdg = len(self.readings)
        for k in range(num_rdg):
            num_dipole = len(self.readings[k].Vdp)
            for j in range(num_dipole):
                if reject is "app_rho":
                    if self.readings[k].Vdp[j].flagRho == "Accept":
                        error_list.append(self.readings[k].Vdp[j].Mx_err)
                elif reject is "app_mx":
                    if self.readings[k].Vdp[j].flagMx == "Accept":
                        error_list.append(self.readings[k].Vdp[j].Mx_err)
                else:
                    error_list.append(self.readings[k].Vdp[j].Mx_err)

        return np.vstack(error_list)

    def getVpDivInErrors(self, reject=None):
        """
        Exports all the eta data

        Output:
        numpy array [value]

        """
        error_list = []
        num_rdg = len(self.readings)
        for k in range(num_rdg):
            num_dipole = len(self.readings[k].Vdp)
            for j in range(num_dipole):
                if reject == 'app_mx':
                    if self.readings[k].Vdp[j].flagMx == "Accept":
                        err = (self.readings[k].Vdp[j].Vp_err + self.readings[k].Vdp[j].In_err) / 100.
                        error_list.append(err)
                elif reject == 'app_rho':
                    if self.readings[k].Vdp[j].flagRho == "Accept":
                        err = (self.readings[k].Vdp[j].Vp_err + self.readings[k].Vdp[j].In_err) / 100.
                        error_list.append(err)
                else:
                    err = (self.readings[k].Vdp[j].Vp_err + self.readings[k].Vdp[j].In_err) / 100.
                    error_list.append(err)

        return np.vstack(error_list)

    def getActiveIndicies(self, reject=None):
        """
        exports an array containing index of each available dipole

        Output:
        list []
        """
        index_list = []
        num_rdg = len(self.readings)
        if reject is not None:
            for k in range(num_rdg):
                num_dipole = len(self.readings[k].Vdp)
                for j in range(num_dipole):
                    if reject is "app_rho":
                        if self.readings[k].Vdp[j].flagRho == "Accept":
                            index = np.array([k, j])
                            index_list.append(index)
                    elif reject is "app_mx":
                        if self.readings[k].Vdp[j].flagMx == "Accept":
                            index = np.array([k, j])
                            index_list.append(index)
        else:
            print("Please choose which indicies you are looking for")

        return np.vstack(index_list)

    def getGeometricFactor(self, reject=None):
        """
        Exports all the geometry factor data

        Output:
        numpy array [value]

        """
        k_list = []
        num_rdg = len(self.readings)
        for k in range(num_rdg):
            num_dipole = len(self.readings[k].Vdp)
            for j in range(num_dipole):
                if reject == 'app_rho':
                    if (self.readings[k].Vdp[j].flagRho == "Accept"):
                        # k_a = self.readings[k].Vdp[j].K
                        k_a = self.readings[k].Vdp[j].calcGeoFactor(self.readings[k].Idp)
                        try:
                            k_list.append(k_a)
                        except ZeroDivisionError:
                            k_list.append(0);
                elif reject == 'app_mx':
                    if (self.readings[k].Vdp[j].flagMx == "Accept"):
                        k_a = self.readings[k].Vdp[j].calcGeoFactor(self.readings[k].Idp)
                        try:
                            k_list.append(k_a)
                        except ZeroDivisionError:
                            k_list.append(0);
                else:
                    k_a = self.readings[k].Vdp[j].calcGeoFactor(self.readings[k].Idp)
                    try:
                        k_list.append(k_a)
                    except ZeroDivisionError:
                        k_list.append(0);
        return np.asarray(k_list)

    def getVoltages(self, reject=None):
        """
        Exports all the geometry factor data

        Output:
        numpy array [value]

        """
        k_list = []
        num_rdg = len(self.readings)
        for k in range(num_rdg):
            num_dipole = len(self.readings[k].Vdp)
            for j in range(num_dipole):
                if reject == 'app_mx':
                    if self.readings[k].Vdp[j].flagMx == "Accept":
                        Vp = self.readings[k].Vdp[j].Vp
                        k_a = Vp
                        k_list.append(k_a)
                elif reject == 'app_rho':
                    if self.readings[k].Vdp[j].flagRho == "Accept":
                        Vp = self.readings[k].Vdp[j].Vp
                        k_a = Vp
                        k_list.append(k_a)
                else:
                    Vp = self.readings[k].Vdp[j].Vp
                    k_a = Vp
                    k_list.append(k_a)

        return np.asarray(k_list)

    def reCalcApparentChageability(self, start, stop, reject=None):
        """
        Exports all the apparent chargeability data

        Output:
        numpy array [value]

        """
        chargeability_list = []
        num_rdg = len(self.readings)
        cnt = 0
        window_widths = self.window_width
        for k in range(num_rdg):
            num_dipole = len(self.readings[k].Vdp)
            for j in range(num_dipole):
                if reject == 'app_mx':
                    if self.readings[k].Vdp[j].flagMx == "Accept":
                        mx_a = self.readings[k].Vdp[j].calcMx(window_widths, start, stop)
                        chargeability_list.append(mx_a)
                elif reject == 'app_rho':
                    if self.readings[k].Vdp[j].flagRho == "Accept":
                        mx_a = self.readings[k].Vdp[j].calcMx(window_widths, start, stop)
                        chargeability_list.append(mx_a)
                else:
                    mx_a = self.readings[k].Vdp[j].calcMx(window_widths, start, stop)
                    chargeability_list.append(mx_a)

        return np.asarray(chargeability_list)

    def getApparentChageability(self, reject=None):
        """
        Exports all the apparent chargeability data

        Output:
        numpy array [value]

        """
        chargeability_list = []
        num_rdg = len(self.readings)
        cnt = 0
        for k in range(num_rdg):
            num_dipole = len(self.readings[k].Vdp)
            for j in range(num_dipole):
                if reject == 'app_mx':
                    if self.readings[k].Vdp[j].flagMx == "Accept":
                        vs = np.sum(self.readings[k].Vdp[j].Vs)
                        Vs = self.readings[k].Vdp[j].Vs
                        Vp = self.readings[k].Vdp[j].Vp
                        mx_a = (self.readings[k].Vdp[j].Mx * (self.readings[k].Vdp[j].Vp / 1e3))
                        Vp = self.readings[k].Vdp[j].Vp
                        # Vs = self.readings[k].Vdp[j].Vs
                        cnt = cnt + 1
                        mx_a = (mx_a) / (Vp / 1e3)
                            # print("{0} VS: {1} & VP: {2}".format(cnt, vs, Vp))
                        chargeability_list.append(mx_a)
                elif reject == 'app_rho':
                    if self.readings[k].Vdp[j].flagRho == "Accept":
                        vs = np.sum(self.readings[k].Vdp[j].Vs)
                        Vs = self.readings[k].Vdp[j].Vs
                        Vp = self.readings[k].Vdp[j].Vp
                        mx_a = (self.readings[k].Vdp[j].Mx * (self.readings[k].Vdp[j].Vp / 1e3))
                        Vp = self.readings[k].Vdp[j].Vp
                        # Vs = self.readings[k].Vdp[j].Vs
                        cnt = cnt + 1
                        mx_a = (mx_a) / (Vp / 1e3)
                            # print("{0} VS: {1} & VP: {2}".format(cnt, vs, Vp))
                        chargeability_list.append(mx_a)
                else:
                    vs = np.sum(self.readings[k].Vdp[j].Vs)
                    Vs = self.readings[k].Vdp[j].Vs
                    Vp = self.readings[k].Vdp[j].Vp
                    mx_a = (self.readings[k].Vdp[j].Mx * (self.readings[k].Vdp[j].Vp / 1e3))
                    Vp = self.readings[k].Vdp[j].Vp
                    # Vs = self.readings[k].Vdp[j].Vs
                    cnt = cnt + 1
                    mx_a = (mx_a) / (Vp / 1e3)
                        # print("{0} VS: {1} & VP: {2}".format(cnt, vs, Vp))
                    chargeability_list.append(mx_a)

        return np.asarray(chargeability_list)

    def getNullFactors(self, reject=None):
        """
        Exports all the coupling factor data

        Output:
        numpy array [value]

        """
        null_list = []
        num_rdg = len(self.readings)
        for rdg in range(num_rdg):
            num_dipole = len(self.readings[rdg].Vdp)
            for dp in range(num_dipole):
                if reject == 'app_rho':
                    if self.readings[rdg].Vdp[dp].flagRho == "Accept":
                        null_list.append(self.readings[rdg].Vdp[dp].coupling)
                elif reject == 'app_mx':
                    if self.readings[rdg].Vdp[dp].flagMx == "Accept":
                        null_list.append(self.readings[rdg].Vdp[dp].coupling)
                else:
                        null_list.append(self.readings[rdg].Vdp[dp].coupling)

        return np.asarray(null_list)

    def getCurrents(self, reject=None):
        """
        Exports all the In data

        Output:
        numpy array [value]

        """
        in_list = []
        num_rdg = len(self.readings)
        for k in range(num_rdg):
            num_dipole = len(self.readings[k].Vdp)
            for j in range(num_dipole):
                if reject == 'app_rho':
                    if self.readings[k].Vdp[j].flagRho == "Accept":
                        In = np.abs(self.readings[k].Vdp[j].In)
                        in_list.append(In)
                elif reject == 'app_mx':
                    if self.readings[k].Vdp[j].flagMx == "Accept":
                        In = np.abs(self.readings[k].Vdp[j].In)
                        in_list.append(In)
                else:
                    In = np.abs(self.readings[k].Vdp[j].In)
                    in_list.append(In)

        return np.asarray(in_list)

    def getAspacings(self, local=False, direction=None, reject=None):
        
        a_list = []
        num_rdg = len(self.readings)
        if local:
            if direction is not None:
                for k in range(num_rdg):
                    num_dipole = len(self.readings[k].Vdp)
                    for j in range(num_dipole):
                        if reject == 'app_rho':
                            if self.readings[k].Vdp[j].flagRho == "Accept":
                                sep = np.abs(self.readings[k].Vdp[j].getAseperationLocal(direction=direction))
                                a_list.append(sep)
                        if reject == 'app_mx':
                            if self.readings[k].Vdp[j].flagMx == "Accept":
                                sep = np.abs(self.readings[k].Vdp[j].getAseperationLocal(direction=direction))
                                a_list.append(sep)
                        else:
                            sep = np.abs(self.readings[k].Vdp[j].getAseperationLocal(direction=direction))
                            a_list.append(sep)
            else:
                print("[OUTPUT]  Please define line direction for local option!!!!!!!!!")
        else:
            for k in range(num_rdg):
                num_dipole = len(self.readings[k].Vdp)
                for j in range(num_dipole):
                    if reject == 'app_rho':
                        if self.readings[k].Vdp[j].flagRho == "Accept":
                            sep = np.abs(self.readings[k].Vdp[j].getAseperation())
                            a_list.append(sep)
                    elif reject == 'app_mx':
                        if self.readings[k].Vdp[j].flagMx == "Accept":
                            sep = np.abs(self.readings[k].Vdp[j].getAseperation())
                            a_list.append(sep)
                    else:
                        sep = np.abs(self.readings[k].Vdp[j].getAseperation())
                        a_list.append(sep)
        return a_list

    def getSources(self, dipole=False):
        """
        Exports all the tx locations to a numpy array

        Output:
        numpy array [x0,y0,z0,x1,y1,z1] (c1,c2)

        """

        src_list = []
        num_rdg = len(self.readings)
        if not dipole:
            for k in range(num_rdg):
                tx = np.array([self.readings[k].Idp.Tx1East,
                              self.readings[k].Idp.Tx1North,
                              self.readings[k].Idp.Tx1Elev,
                              self.readings[k].Idp.Tx2East,
                              self.readings[k].Idp.Tx2North,
                              self.readings[k].Idp.Tx2Elev])
                src_list.append(tx)
        else:
            for k in range(num_rdg):
                num_dipole = len(self.readings[k].Vdp)
                for j in range(num_dipole):
                    tx = np.array([self.readings[k].Idp.Tx1East,
                                  self.readings[k].Idp.Tx1North,
                                  self.readings[k].Idp.Tx1Elev,
                                  self.readings[k].Idp.Tx2East,
                                  self.readings[k].Idp.Tx2North,
                                  self.readings[k].Idp.Tx2Elev])
                    src_list.append(tx)
        # number_of_src = len(src_list)
        # src = np.zeros((number_of_src, 6))
        # for tx in range(number_of_src):
        #     src[tx, :] = src_list[tx]
        return np.asarray(src_list)

    def getSources2(self, dipole=False):
        """
        Exports all the tx locations to a numpy array

        Output:
        numpy array [x0,y0,z0,x1,y1,z1] (c1,c2)

        """

        tx = []
        num_rdg = len(self.readings)
        if not dipole:
            for k in range(num_rdg):
                tx = np.array([self.readings[k].Idp.Tx1East,
                              self.readings[k].Idp.Tx1North,
                              self.readings[k].Idp.Tx1Elev])
                src_list.append(tx)
        else:
            for k in range(num_rdg):
                num_dipole = len(self.readings[k].Vdp)
                for j in range(num_dipole):
                    # tx = np.array([self.readings[k].Idp.Tx1East,
                    #               self.readings[k].Idp.Tx1North,
                    #               self.readings[k].Idp.Tx1Elev])
                    if self.readings[k].Vdp[j].flagMx == "Accept":
                        tx.append(np.array([self.readings[k].Idp.Tx1East,
                                            self.readings[k].Idp.Tx1North]))
                src_list = np.vstack(tx)
        # number_of_src = len(src_list)
        # src = np.zeros((number_of_src, 6))
        # for tx in range(number_of_src):
        #     src[tx, :] = src_list[tx]
        return np.asarray(src_list)

    def getDipoles(self):
        """
        Exports all the tx locations to a numpy array

        Output:
        numpy array [x0,y0,z0,x1,y1,z1] (p1,p2)

        """

        dipole_list = []
        num_rdg = len(self.readings)
        for k in range(num_rdg):
            num_dipole = len(self.readings[k].Vdp)
            for j in range(num_dipole):
                rx = np.array([self.readings[k].Vdp[j].Rx1East,
                              self.readings[k].Vdp[j].Rx1North,
                              self.readings[k].Vdp[j].Rx1Elev,
                              self.readings[k].Vdp[j].Rx2East,
                              self.readings[k].Vdp[j].Rx2North,
                              self.readings[k].Vdp[j].Rx2Elev])
                dipole_list.append(rx)
        number_of_dipoles = len(dipole_list)
        dipoles = np.zeros((number_of_dipoles, 6))
        for rx in range(number_of_dipoles):
            dipoles[rx, :] = dipole_list[rx]
        return dipoles

    def getNumberOfDataPoints(self, active=False):
        """
        returns number of dipoles in the dataset
        """
        num_rdg = len(self.readings)
        count = 0
        for k in range(num_rdg):
            num_dipole = len(self.readings[k].Vdp)
            for j in range(num_dipole):
                # if active:
                #     if self.readings[k].Vdp[j].flagRho
                count += 1
        return count

    def CompareDatabase2RecalcValues(self, mx_start, mx_end, reject=None):
        calc_mx = self.reCalcApparentChageability(mx_start, mx_end, reject=reject)
        mx = self.getApparentChageability(reject=reject)
        calc_rho = self.getApparentResistivity(reject=reject)
        rho = self.getApparentResistivity(reject=reject, calc=False)
        fig = plt.figure()
        ax1 = fig.add_subplot(221)
        ax2 = fig.add_subplot(222)
        ax3 = fig.add_subplot(223)
        ax4 = fig.add_subplot(224)

        ax1.set_title(r'$\rho$')
        ax2.set_title(r'$\rho$' + ' diff')
        ax3.set_title(r'$\eta$')
        ax4.set_title(r'$\eta$' + ' diff')

        ax1.plot(mx, '.')
        ax1.plot(calc_mx, '.r', markersize=1)
        ax2.plot(calc_mx - mx, 'or')
        ax3.plot(rho, '.')
        ax3.plot(calc_rho, '.r', markersize=1)
        ax4.plot(calc_rho - rho, 'or')
        plt.show()

    def createDcSurvey(self, data_type, ip_type=None):
        """
        Loads a dias data file to a SimPEG "srcList" class

        Input:
        datatype = Choose either IP or DC

        Note: elevation is +ve for simPEG inversion

        """
        doff = 0                                      # in case offset is require
        srcLists = []                                 # Injections + dipoles
        data = []                                     # data from file
        d_weights = []                                # weights for the data
        num_rdg = len(self.readings)
        minE = self.readings[0].Idp.Tx2East
        minN = self.readings[0].Idp.Tx2North
        maxN = minN
        maxE = minE
        for k in range(num_rdg):
            num_dipole = len(self.readings[k].Vdp)
            num_dipole_count = 0
            for i in range(num_dipole):
                if data_type == "DC":
                    if (self.readings[k].Vdp[i].flagRho == "Accept"):
                        num_dipole_count += 1
                if data_type == "IP":
                    if self.readings[k].Vdp[i].flagMx == "Accept":
                        num_dipole_count += 1
            rx = np.zeros((num_dipole_count, 6))
            tx = np.array([self.readings[k].Idp.Tx1East,
                          self.readings[k].Idp.Tx1North,
                          self.readings[k].Idp.Tx1Elev - doff,
                          # 0,
                          self.readings[k].Idp.Tx2East,
                          self.readings[k].Idp.Tx2North,
                          self.readings[k].Idp.Tx2Elev - doff])
                          # 0])
            if self.readings[k].Idp.Tx1East > maxE:
                maxE = self.readings[k].Idp.Tx1East
            if self.readings[k].Idp.Tx1East < minE:
                minE = self.readings[k].Idp.Tx1East
            if self.readings[k].Idp.Tx1North > maxN:
                maxN = self.readings[k].Idp.Tx1North
            if self.readings[k].Idp.Tx1North < minN:
                minN = self.readings[k].Idp.Tx1North
            cnt = 0
            for i in range(num_dipole):
                if data_type == "DC":
                    if (self.readings[k].Vdp[i].flagRho == "Accept"):
                        rx[cnt, :] = [self.readings[k].Vdp[i].Rx1East,
                                      self.readings[k].Vdp[i].Rx1North,
                                      self.readings[k].Vdp[i].Rx1Elev - doff,
                                      # 0,
                                      self.readings[k].Vdp[i].Rx2East,
                                      self.readings[k].Vdp[i].Rx2North,
                                      self.readings[k].Vdp[i].Rx2Elev - doff]
                                      # 0]
                        Vp = self.readings[k].Vdp[i].Vp
                        data.append(Vp / self.readings[k].Idp.In)
                        d_weights.append((self.readings[k].Vdp[i].Vp_err +
                                          self.readings[k].Idp.In_err) / 100.)
                        cnt += 1
                if data_type == "IP":
                    if self.readings[k].Vdp[i].flagMx == "Accept":
                        rx[cnt, :] = [self.readings[k].Vdp[i].Rx1East,
                                      self.readings[k].Vdp[i].Rx1North,
                                      self.readings[k].Vdp[i].Rx1Elev - doff,
                                      self.readings[k].Vdp[i].Rx2East,
                                      self.readings[k].Vdp[i].Rx2North,
                                      self.readings[k].Vdp[i].Rx2Elev - doff]
                        if ip_type is None:
                            data.append(self.readings[k].Vdp[i].Mx)
                        elif ip_type == "decay":
                            Vp = self.readings[k].Vdp[i].Vp
                            data.append(self.readings[k].Vdp[i].Vs)
                        else:
                            Vp = self.readings[k].Vdp[i].Vp
                            data.append(self.readings[k].Vdp[i].Mx *
                                        ((Vp / 1e3) / self.readings[k].Idp.In))

                        d_weights.append((self.readings[k].Vdp[i].Mx *
                                          (self.readings[k].Vdp[i].Mx_err /
                                           100.0)))
                        cnt += 1

            Rx = DC.Rx.Dipole(rx[:, :3], rx[:, 3:])    # create dipole list
            # srcLists.append(DC.Src.Pole([Rx], tx[:3]))
            srcLists.append(DC.Src.Dipole([Rx], tx[:3], tx[3:]))

        survey = DC.SurveyDC.Survey(srcLists)          # creates the survey
        survey.dobs = np.float64(np.asarray(data))                 # assigns data
        survey.std = np.asarray(d_weights)             # assign data weights
        survey.eps = 0.001

        return survey


    def plotGpsOverDatabaseLocations(self, gps_input=None):
        if gps_input is not None:
            tx = self.getSources()
            rx = self.getDipoles()
            plt.plot(tx[:, 0], tx[:, 1], 'oc')
            plt.plot(tx[:, 3], tx[:, 4], 'oc')
            plt.plot(rx[:, 0], rx[:, 1], 'oc')
            plt.plot(gps_input[0, :], gps_input[1, :], '+k')
            plt.xlabel('Easting (m)')
            plt.ylabel('Northing (m)')
            plt.title('Database Locations & GPS View')
            plt.show()
        else:
            print("[INFO] Please provide gps input data!!!!")


    def checkDipoleAppRhoPolarityPerReading(self, num_readings=None,
                                            gps_locations=None,
                                            dipole_dipole=False):
        if num_readings is not None:
            num_rdg = num_readings
            for rdg in range(int(num_rdg)):
                self.readings[rdg].createNodeDB()
                if gps_locations is not None:
                    plt.plot(gps_locations[0, :], gps_locations[1, :], '+k')
                if dipole_dipole:
                    tx1x = self.readings[rdg].Idp.Tx1East
                    tx1y = self.readings[rdg].Idp.Tx1North
                    tx2x = self.readings[rdg].Idp.Tx2East
                    tx2y = self.readings[rdg].Idp.Tx2North
                    plt.plot([tx1x, tx2x], [tx1y, tx2y], '-dk')
                else:
                    tx1x = self.readings[rdg].Idp.Tx1East
                    tx1y = self.readings[rdg].Idp.Tx1North
                    plt.plot(tx1x, tx1y, 'dk')

                for idx in range(len(self.readings[rdg].node_db)):
                    x = self.readings[rdg].node_locs[idx][0]
                    y = self.readings[rdg].node_locs[idx][1]
                    plt.plot(x, y, 'k*')
                    plt.text(x, y, '%s' % self.readings[rdg].node_db[idx])
                plt.title("rec: {0}".format(self.readings[rdg].MemNumber))
                for dp in range(len(self.readings[rdg].Vdp)):
                    rx1x = self.readings[rdg].Vdp[dp].Rx1East
                    rx1y = self.readings[rdg].Vdp[dp].Rx1North
                    rx2x = self.readings[rdg].Vdp[dp].Rx2East
                    rx2y = self.readings[rdg].Vdp[dp].Rx2North
                    In = self.readings[rdg].Idp
                    if self.readings[rdg].Vdp[dp].calcRho(In) > 0:
                        plt.plot([rx1x, rx2x], [rx1y, rx2y], '-ob')
                    else:
                        plt.plot([rx1x, rx2x], [rx1y, rx2y], '-or')
                    
                plt.show()
        else:
            num_rdg = len(self.readings)
            for rdg in range(int(num_rdg)):
                self.readings[rdg].createNodeDB()
                if gps_locations is not None:
                    plt.plot(gps_locations[0, :], gps_locations[1, :], '+k')
                for idx in range(len(self.readings[rdg].node_db)):
                    x = self.readings[rdg].node_locs[idx][0]
                    y = self.readings[rdg].node_locs[idx][1]
                    plt.plot(x, y, 'k*')
                    plt.text(x, y, '%s' % self.readings[rdg].node_db[idx])
                plt.title("rec: {0}".format(self.readings[rdg].MemNumber))
                for dp in range(len(self.readings[rdg].Vdp)):
                    rx1x = self.readings[rdg].Vdp[dp].Rx1East
                    rx1y = self.readings[rdg].Vdp[dp].Rx1North
                    rx2x = self.readings[rdg].Vdp[dp].Rx2East
                    rx2y = self.readings[rdg].Vdp[dp].Rx2North
                    In = self.readings[rdg].Idp
                    if self.readings[rdg].Vdp[dp].calcRho(In) > 0:
                        plt.plot([rx1x, rx2x], [rx1y, rx2y], '-ob')
                    else:
                        plt.plot([rx1x, rx2x], [rx1y, rx2y], '-or')
                    
                plt.show()
    

    def plotHistogramOfDataVitals(self, reject='app_rho',
                                  log_rho=False, log_vp=False):
        fig = plt.figure()
        ax1 = fig.add_subplot(241)
        ax2 = fig.add_subplot(242)
        ax3 = fig.add_subplot(243)
        ax4 = fig.add_subplot(244)
        ax5 = fig.add_subplot(245)
        ax6 = fig.add_subplot(246)
        ax7 = fig.add_subplot(247)
        ax8 = fig.add_subplot(248)

        ax1.set_title(r'$\rho$')
        ax2.set_title("error V/I")
        ax3.set_title("K")
        ax4.set_title("V")
        ax5.set_title(r'$\eta$')
        ax6.set_title("Null")
        ax7.set_title("I")
        ax8.set_title("a")

        rho = self.getApparentResistivity(reject=reject)
        vp_i_error = self.getVpDivInErrors(reject=reject)
        gk = self.getGeometricFactor(reject=reject)
        voltages = self.getVoltages(reject=reject)
        mx = self.getApparentChageability(reject=reject)
        null = self.getNullFactors(reject=reject)
        current = self.getCurrents(reject=reject)
        a = self.getAspacings(reject=reject)

        if log_rho:
            ax1.hist(np.log(np.abs(rho)), 50)
            ax1.set_title("log(" + r'$\rho$' + ")")
        else:
            ax1.hist(rho, 50)
        ax2.hist(vp_i_error, 25)
        ax3.hist(gk, 50)
        if log_vp:
            ax4.hist(np.log(voltages), 50)
            ax4.set_title("log(V)")
        else:
            ax4.hist(voltages, 50)
        ax5.hist(mx, 50)
        ax6.hist(null, 50)
        ax7.hist(current, 50)
        ax8.hist(a)

        plt.show()

class Jreadtxtline:
    """
    Class specifically for reading a line of text from a dias file

    """

    def __init__(self, hdrLine, dataLine):
        # make a structure with the header inputs
        self.Vs = []
        for n in range(len(hdrLine)):
            if hdrLine[n].find("Vs") == 0:
                try:
                    self.Vs.append(float(dataLine[n]))
                except ValueError:
                    self.Vs.append(float(-99.9))

            else:
                setattr(self, hdrLine[n], dataLine[n])
