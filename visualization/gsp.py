# adapted from Jay S. Stanley III, Yale University, Fall 2018
# CBB 555 -- Problem Set 2
#
# This script contains functions for implementing graph clustering and signal processing.
#

import numpy as np
import numpy.matlib
import codecs
import sys
import json
import scipy
from scipy import sparse
from scipy.spatial.distance import cdist, pdist, squareform
from scipy.linalg import eigh
from numpy.linalg import norm

from sklearn.cluster import KMeans

def load_json_files(file_path):
    '''
    Loads data from a json file

    Inputs:
        file_path   the path of the .json file that you want to read in

    Outputs:
        my_array    this is a numpy array if data is numeric, it's a list if it's a string

    '''

    #  load data from json file
    obj_text = codecs.open(file_path, 'r', encoding='utf-8').read()
    b_new = json.loads(obj_text)

    # FIXED BUG: simple (but not elegant) hack to make this work for both python versions
    # Python 2 ==> reads the json as unicode so we have to change it
    # Python 3 ==> is smarter and reads the json as ascii

    if sys.version_info[0] is 2:
        # if a string, then returns list of strings
        if type(b_new[0]) is unicode:
            my_array = [b.encode('ascii', 'ignore')
                        for b in b_new]  # translate from unicode to ascii

        # otherwise, it's assumed to be numeric and returns numpy array
        else:
            my_array = np.array(b_new)

    else:
        # if a string, then returns list of strings
        if type(b_new[0]) is str:
            my_array = b_new  # just return

        # otherwise, it's assumed to be numeric and returns numpy array
        else:
            my_array = np.array(b_new)

    return my_array


def gaussian_kernel(X, kernel_type="gaussian", sigma=3.0, k=5):
    """gaussian_kernel: Build an adjacency matrix for data using a Gaussian kernel
    Args:
        X (N x d np.ndarray): Input data
        kernel_type: "gaussian" or "adaptive". Controls bandwidth
        sigma (float): Scalar kernel bandwidth
        k (integer): nearest neighbor kernel bandwidth
    Returns:
        W (N x N np.ndarray): Weight/adjacency matrix induced from X
    """
    _g = "gaussian"
    _a = "adaptive"

    kernel_type = kernel_type.lower()
    D = squareform(pdist(X))
    if kernel_type == _g:  # gaussian bandwidth checking
        print("fixed bandwidth specified")

        if not all([type(sigma) is float, sigma > 0]):  # [float, positive]
            print("invalid gaussian bandwidth, using sigma = max(min(D)) as bandwidth")
            D_find = D + np.eye(np.size(D, 1)) * 1e15
            sigma = np.max(np.min(D_find, 1))
            del D_find
        sigma = np.ones(np.size(D, 1)) * sigma
    elif kernel_type == _a:  # adaptive bandwidth
        print("adaptive bandwidth specified")

        # [integer, positive, less than the total samples]
        if not all([type(k) is int, k > 0, k < np.size(D, 1)]):
            print("invalid adaptive bandwidth, using k=5 as bandwidth")
            k = 5

        knnDST = np.sort(D, axis=1)  # sorted neighbor distances
        sigma = knnDST[:, k]  # k-nn neighbor. 0 is self.
        del knnDST

    W = ((D**2) / sigma[:, np.newaxis]).T
    W = np.exp(-1 * (W))
    W = (W + W.T) / 2  # symmetrize
    W = W - np.eye(W.shape[0])  # remove the diagonal
    return W


def estimate_lmax(L):
    """estimate_lmax:
     Estimate the maximum eigenvalue of the Laplacian.
     borrowed from pygsp.graphs
    Args:
        L (N x N np.ndarray): Matrix to estimate lmax

    Returns:
        lmax (np.float64): estimate of maximum eigenvalue
    """

    N = L.shape[0]
    lmax = sparse.linalg.eigsh(L, k=1, tol=5e-3,
                               ncv=min(N, 10),
                               return_eigenvectors=False)
    lmax = lmax[0]
    lmax *= 1.01
    return lmax

def compute_jackson_cheby_coeff(filter_bounds, delta_lambda, m):
    # compute polynomial approximation coefficients
    if delta_lambda[0] > filter_bounds[0] or delta_lambda[1] < filter_bounds[1]:
        raise("Bounds of the filter are out of the lambda values")
    elif delta_lambda[0] > delta_lambda[1]:
        raise("lambda_min is greater than lambda_max")

    # Scaling and translating to standard cheby interval
    a1 = (delta_lambda[1]-delta_lambda[0])/2
    a2 = (delta_lambda[1]+delta_lambda[0])/2

    # Scaling bounds of the band pass according to lrange
    filter_bounds[0] = (filter_bounds[0]-a2)/a1
    filter_bounds[1] = (filter_bounds[1]-a2)/a1
    # First compute cheby coeffs
    ch = np.arange(float(m+1))
    ch[0] = (2/(np.pi))*(np.arccos(filter_bounds[0])-np.arccos(filter_bounds[1]))
    for i in ch[1:]:
        ch[int(i)] = (2/(np.pi * i)) * \
            (np.sin(i * np.arccos(filter_bounds[0])) - np.sin(i * np.arccos(filter_bounds[1])))

    # Then compute jackson coeffs
    jch = np.arange(float(m+1))
    alpha = (np.pi/(m+2))
    for i in jch:
        jch[int(i)] = (1/np.sin(alpha)) * \
            ((1 - i/(m+2)) * np.sin(alpha) * np.cos(i * alpha) +
             (1/(m+2)) * np.cos(alpha) * np.sin(i * alpha))

    # Combine jackson and cheby coeffs
    jch = ch * jch

    return ch, jch

def estimate_lk(L, k, lmax, Nest, epsilon, order):
    """estimate_lk: Estimate the k-th eigenvalue using dichotomous filtering
    Adapted from Tremblay et. al Matlab CSC Toolbox
    http://cscbox.gforge.inria.fr/
    Args:
        L (N x N np.ndarray): Graph Laplacian
        k (integer): Eigenvalue to estimate
        lmax (float): Estimate of maximum eigenvalue
        Nest (integer): Number of estimations to perform
        epsilon (float): Precision
        order (integer): Number of Jackson-chebyshev polynomial coefficients to compute

    Returns:
        lk (float): k-th eigenvalue
    """
    N = L.shape[0]
    lambda_k_est = np.zeros((Nest, 1))
    ns = int(2 * np.round(np.log(N)))
    # Perform Nest on different of set feature vectors
    for ind_est in range(Nest):
        # Random signals (fixed for one estimation)
        s = np.random.randn(N, ns) * 1 / np.sqrt(ns)

        # Search by dichotomy
        counts = 0
        lambda_min = 0
        lambda_max = lmax
        while (counts != k) and ((lambda_max - lambda_min) / lambda_max > epsilon):
            # Middle of the interval
            lambda_mid = (lambda_min + lambda_max) / 2
            # Approximate ideal filter
            ch, jch = compute_jackson_cheby_coeff(
                [0, lambda_mid], [0, lmax], order)
            # Filtering and counting
            X = cheby_op(L, lmax, jch, s)
            counts = np.round(np.sum(X**2))

            if counts > k:
                lambda_max = lambda_mid
            elif counts < k:
                lambda_min = lambda_mid

        # Store result
        lambda_k_est[ind_est] = (lambda_min + lambda_max) / 2

    # Final estimation
    lk = np.mean(lambda_k_est)
    return lk


def cheby_op(L, lmax, c, signal):
    """cheby_op: Chebyshev polynomial of graph Laplacian applied to input signal.
    Used for approximating graph filters. 
    Adapted from pygsp. 
    Args:
        L (N x N np.ndarray): Graph Laplacian
        lmax (float): Estimated maximum eigenvalue of the graph Laplacian.
        c (integer): Chebyshev coefficients
        signal (np.ndarray): Input signal to filter.
    Returns:
        fs (np.ndarray): Filtered signal output

    """
    N = L.shape[0]
    c = np.atleast_2d(c)
    Nscales, M = c.shape

    if M < 2:
        raise TypeError("The coefficients have an invalid shape")

    try:
        Nv = np.shape(signal)[1]
        r = np.zeros((N * Nscales, Nv))
    except IndexError:
        r = np.zeros((N * Nscales))

    a_arange = [0, lmax]

    a1 = float(a_arange[1] - a_arange[0]) / 2.
    a2 = float(a_arange[1] + a_arange[0]) / 2.

    twf_old = signal
    twf_cur = (L.dot(signal) - a2 * signal) / a1

    tmpN = np.arange(N, dtype=int)
    for i in range(Nscales):
        r[tmpN + N * i] = 0.5 * c[i, 0] * twf_old + c[i, 1] * twf_cur

    factor = 2 / a1 * (L - a2 * sparse.eye(N))
    for k in range(2, M):
        twf_new = factor.dot(twf_cur) - twf_old
        for i in range(Nscales):
            r[tmpN + N * i] += c[i, k] * twf_new

        twf_old = twf_cur
        twf_cur = twf_new

    return r


def upsample(s, m, gamma, Hfilt):
    """upsample a signal given an ideal high pass filter
    Adapted from Tremblay et. al Matlab CSC Toolbox
    http://cscbox.gforge.inria.fr/
    Args:
        s (n_m x 1 np.ndarray): Signal to upsample
        m (n_m x 1 np.ndarray): Known indices
        gamma (float): regularization parameter
        Hfilt (N x N): Ideal High pass filterbank matrix 

    Returns:
        s_gt: upsampled signal
    """
    N = Hfilt.shape[0]
    b = sparse.coo_matrix((s, (m, np.zeros((m.size,)))), shape=(N, 1))
    MtM = sparse.coo_matrix((np.ones((m.size,)), (m, m)), shape=(N, N))

    toinv = MtM + gamma * Hfilt
    s_gt = sparse.linalg.cgs(toinv, b.toarray(), maxiter=100)

    return s_gt


# BEGIN PS2 FUNCTIONS

# helper function
def BlockDiagMat(k,N,p,diag=True):
    N_per_k = int(N//k)
    ones_mat = np.ones((N_per_k,N_per_k))
    
    blocks = np.zeros((k,N_per_k,N_per_k))
    
    for i in range(k):
        blocks[i] = ones_mat
        
    if diag == True:
        block_mat = scipy.linalg.block_diag(*list(blocks))
        block_mat = p*block_mat
        
    if diag == False:
        block_mat = scipy.linalg.block_diag(*list(blocks))
        block_mat = np.ones((N,N)) - block_mat
        block_mat = p*block_mat
    
    return block_mat



def sbm(N, k, pij, pii, sigma):
    """sbm: Construct a stochastic block model

    Args:
        N (integer): Graph size
        k (integer): Number of clusters
        pij (float): Probability of intercluster edges
        pii (float): probability of intracluster edges

    Returns:
        A (numpy.array): Adjacency Matrix
        gt (numpy.array): Ground truth cluster labels
        coords(numpy.array): plotting coordinates for the sbm
    """
    
    #CREATE Adjacency matrix
    block_mat = BlockDiagMat(k,N,p=pii,diag=True) + BlockDiagMat(k,N,p=pij,diag=False)

    for i in range(len(block_mat.reshape(-1))):
        block_mat = block_mat.reshape(-1)
        p_for_adj = np.random.uniform()
        if p_for_adj <= block_mat[i]:
            block_mat[i] = 1
        elif p_for_adj > block_mat[i]:
            block_mat[i] = 0

    A = block_mat.reshape(N,N)
    
    #ground truth vector
    N_per_k = N // k
    gt = np.repeat(np.arange(k),N_per_k)
    
    #coords
    # GROUND TRUTH
    N_per_k = int(N//k)
    deg_incr = 360//k
    deg_means = np.zeros(k)
    deg = float(0)

    for ind, inc in enumerate(deg_means):
        deg_means[ind] = deg
        deg += deg_incr

    x,y = np.cos(np.deg2rad(deg_means)), np.sin(np.deg2rad(deg_means))
    mean_coords = np.column_stack(tup=(x,y))

    #print(mean_coords)

    coords = np.zeros((k,N_per_k,2))

    for i in range(k):
        coords[i] = np.random.normal(loc=mean_coords[i,:],
                                  scale=sigma,
                                  size=(N_per_k,2))

    coords = coords.reshape(N,2)
    return A, gt, coords


def L(A, normalized=True):
    """L: compute a graph laplacian

    Args:
        A (N x N np.ndarray): Adjacency matrix of graph
        normalized (bool, optional): Normalized or combinatorial Laplacian

    Returns:
        L (N x N np.ndarray): graph Laplacian
    """


    # degree matrix
    D =  np.diagflat(np.matrix(A).sum(axis=1).T)
    print(D)
    D_neg_1_2 = np.diagflat(np.power(D.diagonal(),-0.5))
    
    #case 1: normalized Laplacian
    if normalized == True:
        L = np.matmul(D_neg_1_2, np.matmul(A,D_neg_1_2))
        L = np.eye(A.shape[0]) - L
    
    #case 2: unnormalized Laplacian
    elif normalized == False:
        L = D - A


    return L


def compute_fourier_basis(L):
    """compute_fourier_basis: Laplacian Diagonalization

    Args:
        L (N x N np.ndarray): graph Laplacian

    Returns:
        e (N x 1 np.ndarray): graph Laplacian eigenvalues
        psi (N x N np.ndarray): graph Laplacian eigenvectors
    """
    
    #compute eigenvalues and eigenvectors
    e,psi = np.linalg.eigh(L)
    
    
    return e, psi


def gft(s, psi):
    """gft: Graph Fourier Transform (GFT)

    Args:
        s (N x d np.ndarray): Matrix of graph signals.  Each column is a signal.
        psi (N x N np.ndarray): graph Laplacian eigenvectors
    Returns:
        shat (N x d np.ndarray): GFT of the data
    """
    
    # graph signals = feature matrix = s
    
    # convert to frequency domain
    shat = np.matmul(psi.T,s)
    
    return shat


def filterbank_matrix(psi, e, h, c=1):
    """filterbank_matrix: build a filter matrix using the input filter h

    Args:
        psi (N x N np.ndarray): graph Laplacian eigenvectors
        e (N x 1 np.ndarray): graph Laplacian eigenvalues
        h (function handle): A function that takes in eigenvalues
        and returns values in the interval (0,1)

    Returns:
        H (N x N np.ndarray): Filter matrix that can be used in the form
        filtered_s = H@s
    """
    
    #choose filter
    sig = 'sigmoid'
    cos = 'cosine'
    sin = 'sine'
    cut = 'cutoff'
    neg_exp = 'neg_exp'
    
    if h == sig:
        print("h function chosen as sigmoid")
        h_fun = lambda t: 1/(1+np.exp(-t))
        
    elif h == cos:
        print("h function chosen as cosine")
        h_fun = lambda t: np.cos(t)
        
    elif h == sin:
        print("h function chosen as sine")
        h_fun = lambda t: np.sin(t)
        
    elif h == cut:
        print("h function chosen as cutoff")
        print("cutoff chosen as %d" %(c))
        print("setting values above %d to zero" %(c))
        #e[e > c] = 0
        for i in np.arange(np.size(e)):
            if e[i] <= c:
                e[i] = 1
            if e[i] > c:
                e[i] = 0
        h_fun = lambda t: t
        
    elif h == neg_exp:
        print("h function chosen as 1 - exp^-x")
        h_fun = lambda t: 1 - np.exp(-t)
        
    h_fun = np.vectorize(h_fun)
    H = np.matmul(psi, np.matmul(np.diagflat(h_fun(e)), psi.T))
    return H


def kmeans(X, k, nrep=5, itermax=300):
    """kmeans: cluster data into k partitions

    Args:
        X (n x d np.ndarray): input data, rows = points, cols = dimensions
        k (int): Number of clusters to partition
        nrep (int): Number of repetitions to average for final clustering 
        itermax (int): Number of iterations to perform before terminating
    Returns:
        labels (n x 1 np.ndarray): Cluster labels assigned by kmeans
    """
    
    final_centroids = np.zeros((k, X.shape[1]))
    final_labels = np.zeros((X.shape[0],1))
    
    #initalize the distance to beat with each rep
    lowest_total_avg_sum = 1000000
    
    for i in range(nrep):
        #print("\n on {} of {} repetitions".format(i+1,nrep))
        init = kmeans_plusplus(X, k)  # find your initial centroids
        k_clust_members = np.zeros((X.shape[0],k))
        
        
        init_total_avg_sum = 1000000
        current_centroids = init
        

        
        current_iter = 1
        
        while current_iter < itermax:
            
            new_centroids = np.zeros((k, X.shape[1]))
            # step 1: compute distances to centroids - reuturns a n x k matrix
            dist = scipy.spatial.distance.cdist(X,current_centroids)

            # step 2: assign cluster labels based on nearest centroid
            dist_labeled = np.column_stack((dist, (np.argmin(dist,axis=1)+1)))

            # step 3: calculate new centriods based off new members
            for i in range(0, X.shape[0]):
                clust_id = dist_labeled[i,-1]
                k_clust_members[i, int(clust_id)-1] = clust_id

            # loop over all clusters
            k_total_avg_sum = 0
            for i in range(0,k):
                #get current cluster member indices 
                k_clust_ids = np.argwhere(k_clust_members[:,i])
                
                clust_mem_coords = X[k_clust_ids].reshape(-1,X.shape[1])

                #find new centroid coordinates
                sum_x = np.sum(clust_mem_coords[:, 0])
                sum_y = np.sum(clust_mem_coords[:, 1])
                sum_z = 0
                
                #3D case
                if X.shape[1] == 3:  
                    sum_z = np.sum(clust_mem_coords[:, 2])

                #total sum for finding lowest dist centroid
                total_avg_sum = (sum_x + sum_y + sum_z) / X.shape[0]
                k_total_avg_sum += total_avg_sum
                
                #save to current centroid list 
                new_centroids[i,0] = sum_x/clust_mem_coords.shape[0]
                new_centroids[i,1] = sum_y/clust_mem_coords.shape[0]
                
                #3D case
                if X.shape[1] == 3:
                    new_centroids[i,2] = sum_z/clust_mem_coords.shape[0]

   
            centroid_diff = new_centroids - current_centroids
            
            if np.sum(np.abs(centroid_diff)) <= 0.00000:
                #print("convergence found \t iter: {}".format(current_iter))
                #print(new_centroids)
                break

            else: 
                #print('convergence not yet found \t iter: {} \t diff = {}'.format(current_iter,np.sum(centroid_diff)))
                #print(new_centroids)
                current_centroids = new_centroids

            current_iter += 1
            #if current_iter >= itermax:
                #print('convergence not found in {} iterations'.format(current_iter))
                
        #check for minimum distance
        if k_total_avg_sum  < lowest_total_avg_sum:
            
            #if the total average sum of distances (sum of k averages) beats the current
            #min, then those centroids are chose as the final set
            final_centroids = new_centroids
            
            #final set of labels
            final_labels = dist_labeled[:,-1]
            
            #new lowest sum to beat
            lowest_total_avg_sum = k_total_avg_sum
    
    
    labels = final_labels
    
    return labels



def kmeans_plusplus(X, k):
    """kmeans_plusplus: initialization algorithm for kmeans
    Args:
        X (n x d np.ndarray): input data, rows = points, cols = dimensions
        k (int): Number of clusters to partition

    Returns:
        centroids (k x d np.ndarray): centroids for initializing k-means
    """
    
    # step 1: choose a point randomly
    rand_index = np.random.randint(0,X.shape[0])
    #print("initial randomly chosen index point: %d" %(rand_index))
    inital_point = X[rand_index,:]
    
    # step 2: compute distance from this point to all others
    D = squareform(pdist(X))
    dist_vect = D[rand_index]
    
    # step 3: normalize by sum of all distance to get pmf
    dist_pmf = dist_vect/sum(dist_vect)

    # step 4: choose k centriods
    choice_inds = np.random.choice(np.arange(0,X.shape[0]),
                               size = k,
                               replace=False,
                               p=dist_pmf)
    #print("centroid indices: " + str(choice_inds))
    centroids = X[choice_inds,:]
    
    return centroids


def SC(L, k, psi=None, nrep=5, itermax=300, sklearn=False):
    """SC: Perform spectral clustering 
            via the Ng method
    Args:
        L (np.ndarray): Normalized graph Laplacian
        k (integer): number of clusters to compute
        nrep (int): Number of repetitions to average for final clustering
        itermax (int): Number of iterations to perform before terminating
        sklearn (boolean): Flag to use sklearn kmeans to test your algorithm
    Returns:
        labels (N x 1 np.array): Learned cluster labels
    """
    if psi is None:
        # compute the first k elements of the Fourier basis
        # use scipy.linalg.eigh
        psi = np.linalg.eigh(L)[1]
        psisub = psi[:, :k]
        
    else:  # just grab the first k eigenvectors
        psisub = psi[:, :k]

    # normalize your eigenvector rows
    l2_norm = norm(psisub, axis=1, ord=2)
    norm_psi = psisub / l2_norm.reshape(psisub.shape[0],1)
    
    #run k means
    if sklearn:
        labels = KMeans(n_clusters=k, n_init=nrep,
                        max_iter=itermax).fit_predict(psisub)
    else:
        labels = kmeans(X=psisub,k=k, itermax=itermax,nrep=nrep)

    return labels


def CSC(L, k, nsamples=None, nsignals=None, nrep=5, order_lk=10, order_filt=10, gamma=1e-3, itermax=300, sklearn=False):
    """Compressive Spectral Clustering
    https://arxiv.org/abs/1602.02018
    Tremblay et al.
    Args:
        L (N x N np.ndarray): Normalized graph Laplacian
        k (integer): Number of clusters
        nsamples (integer): Amount to downsample
        nsignals (integer): Amount of signals to use for approximation
        nrep (int): Number of repetitions to average for final clustering
        order_lk (integer): Order of the chebyshev polynomial to approx. lk
        order_filt (isanteger):  Order of the chebyshev polynomial to approx. the ideal filter
        gamma (float): amount of regularization to apply when upsampling the labeled points
        itermax (int): Number of iterations to perform before terminating kmeans
        sklearn (boolean): Flag to use sklearn kmeans to test your algorithm

    Returns: 
        labels (N x 1 np.array): Learned cluster labels

    """
    N = L.shape[0]
    
    # set some defaults
    if nsignals is None:
        nsignals = np.round(4 * np.log(N)).astype(int)
    if nsamples is None:
        nsamples = np.round(2 * k * np.log(k)).astype(int)
        
    # estimate the maximum eigenvalue of the graph Laplacian
    lmax = estimate_lmax(L)
    
    # obtain an estimate for the k-th eigenvalue of the Laplacian
    lk = estimate_lk(L, k, lmax, 10, 1e-1, order_lk)
    
    # construct an lk ideal filter
    ch, jch = compute_jackson_cheby_coeff([0, lk], [0, lmax], order_filt)
    
    # generate nsignals random signals
    # generate nsignals using the normal distribution with variance 1/nsignals
    R =  np.random.normal(0,scale=(nsignals**2),size=(N,nsignals))

    HR = cheby_op(L, lmax, jch, R)  # filter the signals using our ideal filter
    
    # normalize the rows with the L2 norm
    l2_norm = norm(HR, axis=1, ord=2)
    HR = HR / l2_norm.reshape(HR.shape[0],-1)
    
    M = np.random.choice(N, nsamples)  # randomly sample nsamples

    # downsample HR using M
    hrm =  HR[M,:]
    
    if sklearn:
        labels = KMeans(n_clusters=k, n_init=nrep,
                        max_iter=itermax).fit_predict(hrm)
    else:
        labels = kmeans(X=hrm, k=k,
                        itermax=itermax,
                       nrep=nrep)
    #labels =  # run your algorithm on the downsampled signals hrm
    #defined above using kmeans
    
    
    # the next steps of the algorithm construct a high pass filter and then recover our true labels
    _, jch = compute_jackson_cheby_coeff(
        [lk, lmax], [0, lmax], order_filt)  # get approximation coefficients
    # evaluate the filter over all the diracs
    Hfilt = cheby_op(L, lmax, jch, np.eye(N))

    # write a function to one hot encode your cluster labels
    # One hot encode your cluster labels. Each cluster gets its own indicator vector
    
    #len of unique labels = should be k 
    n_labels = len(np.unique(labels))

    encoding = dict(zip(  list(np.unique(labels)), list(np.arange(1,n_labels+1))  ))
    labelsmat = np.zeros((len(M),n_labels+1))
    
    for i in range(len(M)):
        cur_label = labels[i]
        col_ind = encoding[cur_label]
        labelsmat[i,int(col_ind)] = 1
    
    #labelsmat =  
    
    outputlabels = np.zeros((N, k))
    
    for i in range(k):
        # upsample each cluster indicator vector
        outputlabels[:, i], info = upsample(labelsmat[:, i], M, gamma, Hfilt)

    # take the cluster with the highest value for each row
    labels = np.zeros(N)
    for i in range(N):
        labels[i] = np.argmax(outputlabels[i,:])

    return labels
