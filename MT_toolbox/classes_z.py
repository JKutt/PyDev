import numpy as np
import pylab as plt
import spectrum_lib
import sys
import glob
import io_lib
import plot_lib

class Z:

    def __init__(self,
                easting=None,
                northing=None,
                elevation=None,
                orientation=None,
                frequencies=None,
                Z=None,
                error=None,
                res_error=None,
                phase_error=None):
        self.easting = easting
        self.northing = northing
        self.elevation = elevation
        self.orientation = orientation
        self.frequencies = frequencies
        self.Z = Z
        self.error = error
        self.res_error = res_error
        self.phase_error = phase_error
        self.mu_0 = 4 * np.pi * 10 ** -1

    def plot_tensor(self):
        plot_lib.plot_tensor(self.frequencies, self.Z, self.error)

    def plot_tensor_res_phase(self):
        self.get_res_phase_error()
        plot_lib.plot_tensor_res_phase(self.frequencies, self.Z, self.res_error, self.phase_error)

    def get_res_phase_error(self):
        self.res_error = np.real(np.asarray([[[2 * self.mu_0 * np.abs(self.Z[i][f][j]) * (self.error[i][f][j]) / self.frequencies[f] / 2. / np.pi
                                                                 for j in range(2)]
                                                                 for f in range(len(self.frequencies))]
                                                                 for i in range(len(self.Z))]))

        self.phase_error = np.real(np.asarray([[[np.arcsin(self.error[i][f][j] / np.abs(self.Z[i][f][j]))
                                                                 for j in range(2)]
                                                                 for f in range(len(self.frequencies))]
                                                                 for i in range(len(self.Z))]))

        for i in range(len(self.Z)):
            for f in range(len(self.frequencies)):
                for j in range(2):
                    if np.isnan(self.phase_error[i][f][j]):
                        self.phase_error[i][f][j] = np.pi / 2.

    def rotate(self):
        print("To do.")

    def get_phase_tensor(self):
        print("To do.")

    def get_invariants(self):
        print("To do.")
