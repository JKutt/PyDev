import numpy as np
import spectrum
import pylab as plt
import scipy.linalg as slin
import time

np.warnings.filterwarnings('ignore')

def get_frequencies(sample_freq, param):
    """

    # Write the frequency array that will be output depending on input parameters.
    # Input:
    # sample_freq: signal sampling frequency
    # param: class of input parameters
    # Output:
    # f_array: frequency array of output MT sounding
    """
    f_array = [0. for i in range(param.nb_increment * (param.nb_reductions + 1))]
    for fact in range(param.nb_reductions + 1):
        freq = np.linspace(0, 1, param.nfft / (param.length_reduction ** fact)) * sample_freq
        for inc in range(param.nb_increment):
            f_array[fact * param.nb_increment + inc] = freq[param.index_first_frequency + inc * param.frequency_increment]

    return f_array


def get_taper(nfft, taper, *args):
    """

    # Return taper with desired parameters
    # Input:
    # nfft: Taper length
    # taper: Taper type
    # *args: depending on the taper
    # Output:
    # taper to be applied to window
    """
    if taper == 'slepian':
        [tapers, eigen] = spectrum.dpss(nfft, args[0], 1)

    return tapers[0]


def apply_taper(window, taper):
    """

    # Return tapered window
    # Input:
    # window: Window to be tapered
    # taper: taper to be applied
    # Output:
    # Taper window
    """
    return np.asarray(window) * np.asarray(taper)


def jackknife(output, input, ref, output_level):
    """

    # Return jackknifed error
    # Input:
    # output: output channels in transfer function system (ie. for MT, electric field)
    # input: input channels in transfer function system (ie. for MT, magnetic field)
    # ref: reference channels in transfer function system (ie. for MT, reference magnetic field)
    # output_level: How chatty you want the code.
    # Output:
    # Jackknifed error on tensor
    """

    z_jackknife = [[] for n in range(len(output))]
    for n in range(len(output)):
        print('\t[INFO] Computing jackknife sample ' + str(n + 1) + '/' + str(len(output)))

        # Removing jackknife sample
        jack_output = output[:].tolist()
        jack_input = [input[i][:].tolist() for i in range(len(input))]
        jack_ref = [ref[i][:].tolist() for i in range(len(ref))]

        del jack_output[n]
        for i in range(len(jack_input)):
            del jack_input[i][n]
        for i in range(len(jack_ref)):
            del jack_ref[i][n]

        jack_output = np.asarray(jack_output)
        jack_input = [np.asarray(jack_input[i]) for i in range(len(jack_input))]
        jack_ref = [np.asarray(jack_ref[i]) for i in range(len(jack_ref))]

        if output_level > 1:
            print('\t[INFO] Computing least-squares solution for first iteration. Number of samples: ' + str(int(len(output))))
        weights = np.asarray([1. for sample in jack_output])
        G, GR, d, weights_matrix = organize_matrices(jack_output, jack_input, jack_ref, weights)
        m = invert_tensor(G, GR, d, weights_matrix)
        scale, abs_error, sum_error = error_on_output(m, d, G, weights)
        count = 0 # counter to prevent infinite loop if unbalanced regression
        while count < 10:
            if output_level > 0:
                print('\t[INFO] Robust regression. Iteration ' + str(count))
            sum_error_pre = sum_error
            # Calculting new weights based on distance error in regression plot
            weights = new_weights(abs_error, scale)
            # Calculating new solution
            G, GR, d, weights_matrix = organize_matrices(jack_output, jack_input, ref, weights)
            m = invert_tensor(G, GR, d, weights_matrix)
            # Calculating new scale and new weighted error
            scale, abs_error, sum_error = error_on_output(m, d, G, weights)
            # Comparison to see if convergence
            if (np.abs(sum_error - sum_error_pre) / sum_error_pre < 0.1):
                break
            else:
                count += 1
        if count >= 10:
            print('\t[WARNING] No convergence found.')
        else:
            if output_level > 0:
                print('\t[INFO] Convergence found.')
            z_jackknife[n] = m

    # Error to be computed here.


def robust_regression(output, input, ref, output_level):
    """

    # Return jackknifed error
    # Input:
    # output: output channels in transfer function system (ie. for MT, electric field)
    # input: input channels in transfer function system (ie. for MT, magnetic field)
    # ref: reference channels in transfer function system (ie. for MT, reference magnetic field)
    # output_level: How chatty you want the code.
    # Output:
    # Tensor
    """
    if output_level > 0:
        print('\t[INFO] Computing least-squares solution for first iteration. Number of samples: ' + str(int(len(output))))
    weights = np.asarray([1. for sample in output])
    G, GR, d, weights_matrix = organize_matrices(output, input, ref, weights)
    m = invert_tensor(G, GR, d, weights_matrix)
    scale, abs_error, sum_error = error_on_output(m, d, G, weights)
    count = 0 # counter to prevent infinite loop if unbalanced regression
    while count < 10:
        if output_level > 0:
            print('\t[INFO] Robust regression. Iteration ' + str(count))
        sum_error_pre = sum_error
        # Calculting new weights based on distance error in regression plot
        weights = new_weights(abs_error, scale)
        # Calculating new solution
        G, GR, d, weights_matrix = organize_matrices(output, input, ref, weights)
        m = invert_tensor(G, GR, d, weights_matrix)
        # Calculating new scale and new weighted error
        scale, abs_error, sum_error = error_on_output(m, d, G, weights)
        # Comparison to see if convergence
        if (np.abs(sum_error - sum_error_pre) / sum_error_pre < 0.1):
            break
        else:
            count += 1
    if count >= 10:
        print('\t[WARNING] No convergence found. Returning NaN.')
        return [np.nan + np.complex(0, 1) * np.nan, np.nan + np.complex(0, 1) * np.nan]
    else:
        if output_level > 0:
            print('\t[INFO] Convergence found.')
        return [m[0][0], m[1][0]]


def new_weights(abs_error, scale):
    """

    # Return weights depending on error scale
    # Input:
    # abs_error = dot product of spread in electric field between measured and predicted
    # scale = here, normalized MAD of the error on electric field
    # Output:
    # Weights to be applied to the system
    """
    weights = []
    for sample_error in abs_error:
        if sample_error[0] / scale <= 1.5:
            weights.append(1)
        else:
            weights.append(1.5 / np.abs(sample_error[0] / scale))
    #weights = np.reshape(weights, (len(weights), 1))

    return weights

def error_on_output(m, d, G, weights):
    """

    # Return error and scale information for the robust regression
    # Input:
    # m: model (ie Z) in the Gm = d problem to solve
    # d: data (ie electric field)
    # G: Predictor (ie magnetic field, 2xN matrix)
    # weights: applied weights to the system
    # Output:
    # scale: MAD/0.44845 to normalize it (ie make it consistent with a std)
    # abs_error: d - Gm
    # sum_error: dot product of abs_error with weights.
    """

    error = np.asarray(d - np.matrix(G) * np.matrix(m))

    abs_error = np.abs(np.conjugate(error) * error)
    scale = np.median(np.abs(abs_error - np.median(abs_error))) / 0.44845
    if scale == 0:
        return np.asarray([[np.nan, np.nan]])

    weights = np.reshape(np.asarray(weights), (len(weights), 1))
    abs_error_weighted = np.abs(np.conjugate(error) * np.asarray(weights) * np.asarray(error))
    sum_error = np.dot(abs_error_weighted[:, 0], abs_error_weighted[:, 0])

    return scale, abs_error, sum_error


def organize_matrices(output, input, ref, weights):
    """

    # Return matrices for linear inversion
    # Input:
    # output: output channels in transfer function system (ie. for MT, electric field)
    # input: input channels in transfer function system (ie. for MT, magnetic field)
    # ref: reference channels in transfer function system (ie. for MT, reference magnetic field)
    # weights: applied weights to the system
    # Output:
    # G: 2xN matrix of the local magnetic field
    # GR: 2xN matrix of the remote magnetic field
    # d: 1xN matrix of the electric field
    # weights_matrix: (NxN) diagonal matrix with weights.
    """

    n_segments = len(output)
    G = np.asarray([[input[i][j] for i in range(len(input))] for j in range(n_segments)])
    GR = np.asarray([[ref[i][j] for i in range(len(ref))] for j in range(n_segments)])
    d = np.asarray([[output[j]] for j in range(n_segments)])


    weights_matrix = [[0. for j in range(n_segments)] for i in range(n_segments)]
    for i in range(n_segments):
        weights_matrix[i][i] = weights[i]

    return G, GR, d, weights_matrix

def invert_tensor(G, GR, d, weights_matrix):
    """

    # Return model issued from linear inversion of least-squares problem
    # Input:
    # G: 2xN matrix of the local magnetic field
    # GR: 2xN matrix of the remote magnetic field
    # d: 1xN matrix of the electric field
    # weights_matrix: (NxN) diagonal matrix with weights.
    # Output:
    # m: model in the linear regression (ie Z)
    """

    cm = np.matrix(GR.transpose().conjugate()) * np.matrix(weights_matrix) * np.matrix(G)
    d = np.matrix(GR.transpose().conjugate()) * np.matrix(weights_matrix) * d

    m = np.asarray(slin.solve(cm, d))

    return m
