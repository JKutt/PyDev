import matplotlib.pylab as plt
import numpy as np

def plot_tensor(frequencies, Z):
    """

    # Plot tensor
    # Input:
    # frequencies: Sampled frequencies in the tensor
    # Z: Tensor
    # Output:
    """


    n_output = len(Z)
    fig, ax = plt.subplots(n_output, 2)
    for channel in range(n_output):
        ax[channel][0].set_xscale('log')
        ax[channel][0].set_yscale('log')
        ax[channel][1].set_xscale('log')
        ax[channel][0].set_xlabel('Frequency (Hz)')
        ax[channel][0].set_ylabel('Module (mV/km/nT)')
        ax[channel][1].set_xlabel('Frequency (Hz)')
        ax[channel][1].set_ylabel('Phase (Rad)')
        ax[channel][0].plot(frequencies, [np.abs(Z[channel][f][0]) for f in range(len(frequencies))], 'ro')
        ax[channel][0].plot(frequencies, [np.abs(Z[channel][f][1]) for f in range(len(frequencies))], 'bo')
        ax[channel][1].plot(frequencies, [np.angle(Z[channel][f][0]) for f in range(len(frequencies))], 'ro')
        ax[channel][1].plot(frequencies, [np.angle(Z[channel][f][1]) for f in range(len(frequencies))], 'bo')
    plt.show()
