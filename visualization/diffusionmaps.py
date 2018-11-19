"""
Diffusion maps 

# Skeleton file by Chris Harshaw, Yale University, Fall 2017
# Adapted from Jay Stanley, Yale University, Fall 2018
# CPSC 453 -- Problem Set 1 

"""

import numpy as np


def compute_distances(X):
    '''
    Constructs a distance matrix from data set, assumes Euclidean distance

    Inputs:
        X       a numpy array of size n x p holding the data set (n observations, p features)

    Outputs:
        D       a numpy array of size n x n containing the euclidean distances between points

    '''
    D_dim = X.T.shape[1]

    D = np.zeros([D_dim, D_dim])

    # looping over all columns for i column
    for i in range(0, mat.T.shape[1]):

        # here we establist the starting indices for j (goes from current i to rest of len)
        for j in range(i, mat.T.shape[1]):

            # get the square root of the sum of the squared difference
            dx = np.sqrt(sum((mat.T[:, i] - mat.T[:, j])**2))

            # enter this into a matrix
            D[i][j] = dx

    # transfom D from upper triangular to symmetric
    D = D + D.T

    # return distance matrix
    return D


def compute_affinity_matrix(D, kernel_type, sigma=None, k=None):
    '''
    Construct an affinity matrix from a distance matrix via gaussian kernel.

    Inputs:
        D               a numpy array of size n x n containing the distances between points
        kernel_type     a string, either "gaussian" or "adaptive".
                            If kernel_type = "gaussian", then sigma must be a positive number
                            If kernel_type = "adaptive", then k must be a positive integer
        sigma           the non-adaptive gaussian kernel parameter
        k               the adaptive kernel parameter

    Outputs:
        W       a numpy array of size n x n that is the affinity matrix

    '''

    if str(kernel_type).lower() == 'gaussian':
        if sigma <= 0:
            raise ValueError('sigma must be a positive number')
        else:
            print("computing affinity matrix with gaussian kernel")
            W = np.exp(-((D) / sigma))
            print("finished!")

    if str(kernel_type).lower() == 'adaptive':
        if k < 0 or type(k) != int:
            raise ValueError('k must be a positive integer')
        else:
            print("computing affinity matrix with adaptive kernel")
            i_nn_vec = np.zeros((D.shape[0], 1))
            j_nn_vec = np.zeros((D.shape[0], 1))

            # form i kNN matrix - diagonal matrix of distance to kth NN
            for i in range(D.shape[0]):
                # get kth nearest neighbor for the ith row (i) - row of distances from xi
                kth_i_nn = np.sort(D[i, :])[k - 1]
                i_nn_vec[i] = kth_i_nn**-2

            i_nn_mat = np.diagflat(i_nn_vec)

            # form j kNN matrix - diagonal matrix of 1/distance to kth NN
            for j in range(D.shape[0]):
                # get kth nearest neighbor for the ith row (i) - row of distances from xi
                kth_j_nn = np.sort(D[:, j])[k - 1]

                j_nn_vec[j] = kth_j_nn**-2

            j_nn_mat = np.diagflat(j_nn_vec)

            # since the adaptive kernel is a sum we can split the operation in two:
            # - the i-matrix where each row is normalized by 2 times its kth nn
            # - the j-matrix where each column is normalized by 2 times its kth nn

            # for the i-matrix, we will multiply the (distance matrix)**2 and then take its transpose
            i_matrix = np.matmul(D, (2 * i_nn_mat))
            i_matrix = np.exp(- i_matrix.T)

            # for the j-matrix, we can just do the same operation as above without taking the transpose at the end
            j_matrix = np.matmul(D, (2 * j_nn_mat))
            j_matrix = np.exp(- j_matrix)
            # print(j_matrix)

            # now we sum these matrices and divide by two
            W = (i_matrix + j_matrix) / 2
            print("finished!")

    if str(kernel_type).lower() not in ['gaussian', 'adaptive']:
        print('kernel type not gaussian or adaptive ')
        print('input matrix returned')
        W = D

    # return the affinity matrix
    return W


def diff_map_info(W):
    '''
    Construct the information necessary to easily construct diffusion map for any t

    Inputs:
        W           a numpy array of size n x n containing the affinities between points

    Outputs:

        diff_vec    a numpy array of size n x n-1 containing the n-1 nontrivial eigenvectors of Markov matrix as columns
        diff_eig    a numpy array of size n-1 containing the n-1 nontrivial eigenvalues of Markov matrix

        We assume the convention that the coordinates in the diffusion vectors are in descending order
        according to eigenvalues.
    '''

    # markov normalize the affinity matrix
    row_sums = W.sum(axis=1)

    # form the row sum diagonal
    D = (row_sums**(-1 / 2)) * np.eye(M=W.shape[0], N=W.shape[1])

    # form the symmetric Markov matrix
    M_s = np.matmul(W, D)
    M_s = np.matmul(D, M_s)

    # eigendecomposition
    diff_eig, diff_vec = (np.linalg.eigh(M_s))

    # normalize eigenvectors
    diff_vec_norm = np.matmul(D, diff_vec) / np.linalg.norm(np.matmul(D, diff_vec), axis=0)

    # remove the trivial eigenvalue
    diff_eig = diff_eig[:-1]

    # remove the trivial eigenvector
    diff_vec_norm = diff_vec_norm[:, :-1]

    # return the info for diffusion maps
    return diff_vec_norm, diff_eig


def get_diff_map(diff_vec, diff_eig, t=1):
    '''
    Construct a diffusion map at t from eigenvalues and eigenvectors of Markov matrix

    Inputs:
        diff_vec    a numpy array of size n x n-1 containing the n-1 nontrivial eigenvectors of Markov matrix as columns
        diff_eig    a numpy array of size n-1 containing the n-1 nontrivial eigenvalues of Markov matrix
        t           diffusion time parameter t

    Outputs:
        diff_map    a numpy array of size n x n-1, the diffusion map defined for t
    '''
    # power eigenvalues
    pow_eigen_vals = diff_eig**t

    # matrix multiply eigenvalues with eigenvectors to get diffusion map
    diff_map = diff_vec * pow_eigen_vals

    return diff_map
