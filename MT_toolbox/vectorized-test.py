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

    weights_matrix = np.zeros(n_segments, n_segments)
    np.fill_diagonal(weights_matrix, weights)
    G_ = np.asarray(input)
    G = np.reshape(G_, (G_.shape[1], 2))
    GR_ = np.asarray(ref)
    GR = np.reshape(GR_, (GR_.shape[1], 2))
    d_ = np.asarray(output)
    d = np.reshape(d_, (d_.shape[1], 2))
    cm = np.matmul(GR.transpose().conjugate(), np.matmul(weights_matrix, np.matrix(G)))
    dm = GR.transpose().conjugate() @ weights_matrix @ d

    return G, GR, d, cm, dm, weights_matrix