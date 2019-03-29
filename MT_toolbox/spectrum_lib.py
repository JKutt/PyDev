import numpy as np
import spectrum
import pylab as plt
import scipy.linalg as slin
import time
from scipy.stats import rayleigh, beta
# import libMatrix
# import libInversion
import scipy

np.warnings.filterwarnings('ignore')

def error_on_output_vectorized(G, d, m, w):


    error = np.zeros((len(d), 1))
    error = d - G @ m
    abs_error = abs(error * np.conjugate(error))
    abs_error_weighted = abs(error * w @ np.conjugate(error))
    sum_error = np.dot(np.reshape(abs_error_weighted, (len(error),)), np.conjugate(np.reshape(abs_error_weighted, (len(error),))))
    return abs_error, sum_error

def organize_matrices_vectorized(output, input, ref, weights):
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

    weights_matrix = np.zeros((n_segments, n_segments))
    np.fill_diagonal(weights_matrix, weights)
    G_ = np.asarray(input)
    G = np.reshape(G_, (G_.shape[0], 2))
    GR_ = np.asarray(ref)
    GR = np.reshape(GR_, (GR_.shape[0], 2))
    d_ = np.asarray(output)
    d = np.reshape(d_, (d_.shape[0], 1))
    cm = np.matmul(GR.transpose().conjugate(), np.matmul(weights_matrix, np.matrix(G)))
    dm = GR.transpose().conjugate() @ weights_matrix @ d

    return G, GR, d, cm, dm, weights_matrix


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
        tapers = np.asarray(tapers)[:, 0]
    if taper == 'hamming':
        tapers = spectrum.window.create_window(nfft, 'hamming')
    if taper == 'window_tukey':
        tapers = spectrum.window.create_window(nfft, 'tuckey')
    if taper == 'null':
        tapers = np.asarray([1. for i in range(nfft)])

    return tapers


def apply_taper(window, taper):
    """

    # Return tapered window
    # Input:
    # window: Window to be tapered
    # taper: taper to be applied
    # Output:
    # Taper window

    freq = np.linspace(0, 1, len(window)) * 2048
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    ax1.set_xscale('log')
    ax2.set_xscale('log')
    ax1.plot(freq, np.abs(np.fft.fft(window)))
    ax2.plot(freq, np.abs(np.fft.fft(np.asarray(window) * np.asarray(taper))))
    plt.show()
    #exit()
    """
    return np.asarray(window) * np.asarray(taper)


def robust_regression(output, input, ref, output_level):

    if output_level > 0:
        print('\t[INFO] Starting robust regression. Number of samples: ' + str(int(len(output))))
        s_time = time.time()

    # Initializing variables for robust regression
    z_error = [0., 0.]
    rquantile = rayleigh.ppf(1 - 1. / len(output), loc=0, scale=1)
    bquantile = beta.ppf(1 - 1. / len(output), 2, len(output) - 2)
    weights = np.asarray([1. for sample in output])
    u_mat = np.diag([np.complex(1, 0) for j in range(len(output))])

    # Preparing arrays for Fortran call
    output = np.asarray(output, order='F')
    input = np.asarray(input, order='F').transpose()
    ref = np.asarray(ref, order='F').transpose()
    weights = np.asarray(weights, order='F')
    hat = hat_matrix(input, u_mat)

    # Preparing matrices
    s_time = time.time()
    #G, GR, d, cm, dm = libMatrix.prepare_matrices(output, input, ref, weights, len(output))
    G, GR, d, cm, dm, _ = organize_matrices_vectorized(output, input, ref, weights)
    if output_level > 2:
        print("\t\t\t[INFO_TIME] Elapsed time for matrix organization: " + str(time.time() - s_time) + ' seconds.')

    # Inverting least-squares solutions
    s_time = time.time()
    m = invert_tensor(cm, dm)
    if output_level > 2:
        print("\t\t\t[INFO_TIME] Elapsed time for matrix inversion: " + str(time.time() - s_time) + ' seconds.')

    # Calculting error on output field
    s_time = time.time()
    abs_error, sum_error = error_on_output_vectorized(G, d, m, weights)
    #abs_error, sum_error = libMatrix.error_on_output(np.asarray(m, order='F').transpose(), d, G, weights, len(output))
    scale = np.median(np.abs(abs_error - np.median(abs_error))) / 0.44845
    if output_level > 2:
        print("\t\t\t[INFO_TIME] Elapsed time for error calculation: " + str(time.time() - s_time) + ' seconds.')

    count = 0 # counter to prevent infinite loop if unbalanced regression
    while count < 10: # Robust regression loop
        if output_level > 0:
            print('\t[INFO] Robust regression. Iteration ' + str(count))
        sum_error_pre = sum_error

        # Calculting new weights based on distance error in regression plot
        weights = new_weights(abs_error, scale)                 # HUBER weights
        l_weights = leverage_weights(hat, u_mat, bquantile)     # LEVERAGE weights
        weights = weights * np.asarray(l_weights)               # Weights combination
        u_mat = np.diag(weights)

        # Calculating new solution
        s_time = time.time()
        G, GR, d, cm, dm, _ = organize_matrices_vectorized(output, input, ref, weights)
        #G, GR, d, cm, dm = libMatrix.prepare_matrices(output, input, ref, weights, len(output))
        if output_level > 2:
            print("\t\t\t[INFO_TIME] Elapsed time for matrix organization: " + str(time.time() - s_time) + ' seconds.')

        # Calculating new solutions using new weights
        s_time = time.time()
        m = invert_tensor(cm, dm)
        if output_level > 2:
            print("\t\t\t[INFO_TIME] Elapsed time for matrix inversion: " + str(time.time() - s_time) + ' seconds.')

        # Calculating new scale and new weighted error
        s_time = time.time()
        abs_error, sum_error = error_on_output_vectorized(G, d, m, weights)
        #abs_error, sum_error = libMatrix.error_on_output(np.asarray(m, order='F').transpose(), d, G, weights, len(output))
        scale = np.median(np.abs(abs_error - np.median(abs_error))) / 0.44845
        if output_level > 2:
            print("\t\t\t[INFO_TIME] Elapsed time for error calculation: " + str(time.time() - s_time) + ' seconds.')

        # Comparison to see if convergence
        if (np.abs(sum_error - sum_error_pre) / sum_error_pre < 0.01):
            break
        else:
            count += 1

        # Update hat matrix
        hat = hat_matrix(input, u_mat)

    if count >= 10:
        print('\t[WARNING] No convergence found. Returning NaN.')
        return [np.nan + np.complex(0, 1) * np.nan, np.nan + np.complex(0, 1) * np.nan]
    else:
        if output_level > 0:
            print('\t[INFO] Convergence found.')

    z_regression = [m[0][0], m[1][0]]

    if output_level > 0:
        print('\t[INFO] Starting jackknife procedure.')
    z_jackknife = [[] for n in range(len(output))]

    # Using full matrices to derive delete-ones
    _, _, _, cm_0, d_0, _ = organize_matrices_vectorized(output, input, ref, weights)
    #cm_0, d_0 = libMatrix.prepare_tensor_inversion(output, input, ref, weights, len(output))
    for n in range(len(output)):
        if output_level > 1:
            print('\t[INFO] Computing jackknife sample ' + str(n + 1) + '/' + str(len(output)))

        # Removing jackknife sample, a bit of algebra shows it's just a substraction
        cm = cm_0 - np.asarray([[ref[n][0] * weights[n] * input[n][0], ref[n][0] * weights[n] * input[n][1]],
                                [ref[n][1] * weights[n] * input[n][0], ref[n][1] * weights[n] * input[n][1]]])
        dm = d_0 - np.asarray([[ref[n][0] * weights[n] * output[n]],[ref[1][0] * weights[n] * output[n]]])

        if output_level > 1:
            print('\t[INFO] Computing least-squares solution for jackknife sample: ' + str(n))

        s_time = time.time()
        m = invert_tensor(cm, dm)
        if output_level > 2:
            print("\t\t\t[INFO_TIME] Elapsed time for matrix inversion: " + str(time.time() - s_time) + ' seconds.')
        z_jackknife[n] = m

    _, z_error[0] = jackknife_statistics([z_jackknife[i][0] for i in range(len(z_jackknife))])
    _, z_error[1] = jackknife_statistics([z_jackknife[i][1] for i in range(len(z_jackknife))])

    return z_regression, z_error



def robust_regression_old(output, input, ref, output_level):
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

    s_time = time.time()
    G, GR, d, cm, dm, weights_matrix = organize_matrices(output, input, ref, weights)
    print("\t[INFO] Elapsed time for matrix organization: " + str(time.time() - s_time) + ' seconds.')

    s_time = time.time()
    m = invert_tensor(cm, dm)
    print("\t[INFO] Elapsed time for matrix inversion: " + str(time.time() - s_time) + ' seconds.')

    s_time = time.time()
    scale, abs_error, sum_error = error_on_output(m, d, G, weights)
    print("\t[INFO] Elapsed time for error calculation: " + str(time.time() - s_time) + ' seconds.')

    count = 0 # counter to prevent infinite loop if unbalanced regression
    while count < 10:
        if output_level > 0:
            print('\t[INFO] Robust regression. Iteration ' + str(count))
        sum_error_pre = sum_error
        # Calculting new weights based on distance error in regression plot
        weights = new_weights(abs_error, scale)
        # Calculating new solution
        G, GR, d, cm, dm, weights_matrix = organize_matrices(output, input, ref, weights)
        m = invert_tensor(cm, dm)
        # Calculating new scale and new weighted error
        scale, abs_error, sum_error = error_on_output(m, d, G, weights)
        # Comparison to see if convergence
        if (np.abs(sum_error - sum_error_pre) / sum_error_pre < 0.01):
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


def jackknife_statistics(population):

    mean = 1. / len(population) * np.sum(population)
    # im sure this could be vectorized / tiled 
    var = (len(population) - 1) / len(population) * np.sum([(sample - mean) ** 2 for sample in population])

    return mean, var


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

def leverage_weights(hat, u, bquantile):
    """
    # Leverage weights based on Chave et Thomson 2014, using the hat matrix
    # Input:
    # hat: matrix of predictor
    # u: combination of Huber and Lev
    """
    weights = [np.complex(1, 0) for i in range(len(hat[:, 0]))]
    for i in range(len(hat[:, 0])):
        y = np.trace(u) * hat[i, i] / 2
        weights[i] = np.real(weights[i]) * np.exp(np.exp(-1 * bquantile ** 2)) *\
                                          np.exp(-1 * np.exp(bquantile * (y - bquantile))) * np.complex(1, 0)
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
    weights_matrix = [[0. for j in range(n_segments)] for i in range(n_segments)]
    for i in range(n_segments):
        weights_matrix[i][i] = weights[i]

    G = np.asarray([[input[i][j] for i in range(len(input))] for j in range(n_segments)])
    GR = np.asarray([[ref[i][j] for i in range(len(ref))] for j in range(n_segments)])
    d = np.asarray([output[j] for j in range(n_segments)])
    weights = np.asarray([weights[j] for j in range(n_segments)])

    weights_matrix = np.asarray(weights_matrix, order='F')

    cm = np.matmul(GR.transpose().conjugate(), np.matmul(weights_matrix, np.matrix(G)))
    dm = GR.transpose().conjugate() @ weights_matrix @ d

    return G, GR, d, cm, dm, weights_matrix

def invert_tensor(cm, d):
    """

    # Return model issued from linear inversion of least-squares problem
    # Input:
    # cm: Left side of least squares regression model containing weights and local predictor
    # d: Right side of least squares regression model, containing weights and predicted variable
    # m: model in the linear regression (ie Z)
    """
    m = np.asarray(slin.solve(cm, d))

    return np.reshape(m, (2, 1))


def hat_matrix(input, weights):
    nC = input[:, 0].size
    hat = np.zeros((nC, nC))

    dot_product_H = np.zeros((2, 2))
    dot_product_H[0, 0] = np.sum(input[:, 0].conjugate() * input[:, 0] * weights.diagonal())
    dot_product_H[0, 1] = np.sum(input[:, 0].conjugate() * input[:, 1] * weights.diagonal())
    dot_product_H[1, 0] = np.sum(input[:, 1].conjugate() * input[:, 0] * weights.diagonal())
    dot_product_H[1, 1] = np.sum(input[:, 1].conjugate() * input[:, 1] * weights.diagonal())

    unit = np.asarray([[np.complex(1, 0), np.complex(0, 0)], [np.complex(0, 0), np.complex(1, 0)]])
    # Inverting dot product
    dot_product_H_inv = slin.solve(dot_product_H, unit)

    arg1 = input[:, 0].conjugate() * dot_product_H_inv[0, 0] + input[:, 1].conjugate() * dot_product_H_inv[0, 1]
    arg2 = input[:, 0].conjugate() * dot_product_H_inv[1, 0] + input[:, 1].conjugate() * dot_product_H_inv[1, 1]
    input_tiled_x = np.tile(input[:, 0], (arg1.size, 1))
    input_tiled_y = np.tile(input[:, 1], (arg2.size, 1))
    hat = input_tiled_x * arg1 + input_tiled_y * arg2

    return hat
