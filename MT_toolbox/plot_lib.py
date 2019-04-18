import matplotlib.pylab as plt
import numpy as np

def plot_tensor(frequencies, Z, error):
    """

    # Plot tensor
    # Input:
    # frequencies: Sampled frequencies in the tensor
    # Z: Tensor
    # Output:
    """

    n_output = len(Z)
    fig, ax = plt.subplots(n_output, 2)
    if n_output > 1:
        for channel in range(n_output):
            ax[channel][0].set_xscale('log')
            ax[channel][1].set_xscale('log')
            if n_output > 1 and channel < 2:
                ax[channel][0].set_yscale('log')
            ax[channel][0].set_xlabel('Frequency (Hz)')
            ax[channel][0].set_ylabel('Module (mV/km/nT)')
            ax[channel][1].set_xlabel('Frequency (Hz)')
            ax[channel][1].set_ylabel('Phase (Rad)')
            ax[channel][0].errorbar(frequencies, [np.abs(Z[channel][f][0]) for f in range(len(frequencies))], yerr=[np.abs(error[channel][f][0]) for f in range(len(frequencies))], fmt='ro')
            ax[channel][0].errorbar(frequencies, [np.abs(Z[channel][f][1]) for f in range(len(frequencies))], yerr=[np.abs(error[channel][f][1]) for f in range(len(frequencies))], fmt='bo')
            ax[channel][1].errorbar(frequencies, [np.angle(Z[channel][f][0]) for f in range(len(frequencies))], yerr=[np.angle(error[channel][f][0]) for f in range(len(frequencies))], fmt='ro')
            ax[channel][1].errorbar(frequencies, [np.angle(Z[channel][f][1]) for f in range(len(frequencies))], yerr=[np.angle(error[channel][f][1]) for f in range(len(frequencies))], fmt='bo')
    else:
        ax[0].set_xscale('log')
        ax[1].set_xscale('log')
        ax[0].set_xlabel('Frequency (Hz)')
        ax[0].set_ylabel('Module (mV/km/nT)')
        ax[1].set_xlabel('Frequency (Hz)')
        ax[1].set_ylabel('Phase (Rad)')
        ax[0].errorbar(frequencies, [np.abs(Z[0][f][0]) for f in range(len(frequencies))], yerr=[np.abs(error[0][f][0]) for f in range(len(frequencies))], fmt='ro')
        ax[0].errorbar(frequencies, [np.abs(Z[0][f][1]) for f in range(len(frequencies))], yerr=[np.abs(error[0][f][1]) for f in range(len(frequencies))], fmt='bo')
        ax[1].errorbar(frequencies, [np.angle(Z[0][f][0]) for f in range(len(frequencies))], yerr=[np.angle(error[0][f][0]) for f in range(len(frequencies))], fmt='ro')
        ax[1].errorbar(frequencies, [np.angle(Z[0][f][1]) for f in range(len(frequencies))], yerr=[np.angle(error[0][f][1]) for f in range(len(frequencies))], fmt='bo')
    plt.show()


def plot_tensor_res_phase(frequencies, Z, res_error, phase_error):
    """

    # Plot tensor RES and PHASE: NEED TO DERIVE ERROR
    # Input:
    # frequencies: Sampled frequencies in the tensor
    # Z: Tensor
    # Output:
    # np.abs(z1[n]) ** 2 / 5. / fScale + np.complex(0, 1) * np.arctan2(np.imag(z1[n]), np.real(z1[n]))
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
        ax[channel][0].errorbar(frequencies, [np.abs(Z[channel][f][0]) ** 2 / 5. / frequencies[f] for f in range(len(frequencies))],
                                yerr=[res_error[channel][f][0] for f in range(len(frequencies))], fmt='ro')
        ax[channel][0].errorbar(frequencies, [np.abs(Z[channel][f][1]) ** 2 / 5. / frequencies[f] for f in range(len(frequencies))],
                                yerr=[res_error[channel][f][1] for f in range(len(frequencies))], fmt='bo')
        ax[channel][1].errorbar(frequencies, [np.angle(Z[channel][f][0]) for f in range(len(frequencies))],
                                yerr=[phase_error[channel][f][1] for f in range(len(frequencies))], fmt='ro')
        ax[channel][1].errorbar(frequencies, [np.angle(Z[channel][f][1]) for f in range(len(frequencies))],
                                yerr=[phase_error[channel][f][1] for f in range(len(frequencies))], fmt='bo')
    plt.show()
