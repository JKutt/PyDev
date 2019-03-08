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
                error=None):
        self.easting = easting
        self.northing = northing
        self.elevation = elevation
        self.orientation = orientation
        self.frequencies = frequencies
        self.Z = Z
        self.error = error

    def plot(self):
        plot_lib.plot_tensor(self.frequencies, self.Z)

    def rotate(self):
        print("To do.")

    def get_phase_tensor(self):
        print("To do.")

    def get_tipper(self):
        print("To do.")

    def get_invariants(self):
        print("To do.")
