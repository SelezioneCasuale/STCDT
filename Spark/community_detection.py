"""
Standalone community detection code for clustering a dataset (usually neural time series)
Written by Mathew H. Evans 2015. Based on 'Spike train community detection'. Humphries, 2011.
"""

import matplotlib.pyplot as plt
import seaborn as sns


# from thunder.viz import Colorize
from numpy import arange, array, ceil, convolve, diff, dot, empty, linspace, load, matrix, median, min,  sum, shape, sqrt, where, zeros
from itertools import cycle
# import pylab
from scipy import stats, spatial
# import scipy.stats
# import scipy.cluster
# import time
from pyspark.mllib.clustering import KMeans
# from math import sqrt
from pyspark import SparkContext

from thunder.utils import pack, getdims, matrices
from thunder.factorization import SVD

def compISIs(self):
    """
	Compute ISIs for each spike train
	"""
    x = where(array(self) >= 1)

    return diff(x)


def convolve_series(self, timescale):
    """
	Convolve each spike train with a Gaussian of width 10*'sigma'.
	Returns the convolved spike train
	"""
    # import scipy.stats

    # Generate the filter and output arrays

    sigma = ceil(timescale).astype(int)

    window = arange(-5*sigma, 5*sigma)

    gaussian = stats.norm.pdf(window, scale=sigma)

    # Make a continuous signal for each train by convolution

    conv_train = convolve(array(self), gaussian)

    return conv_train


def compute_similarity(self, SC):
    """
	Computes pairwise similarity between all spike trains using scipy cosine similarity
	"""
    # import scipy.spatial

    x = self  # np.transpose(self)

    # Find similarity between continuous spike-train vectors
    A = empty([SC.shape[0]])

    for i in range(SC.shape[0]):
        A[i] = spatial.distance.pdist([x, SC[i]], 'cosine')

    # Similarity = 1 - distance
    A = 1 - A

    return A


def zero_diag(k, v):
    """
	Remove diagonal elements of a matrix
	"""
    v[k - 1] = 0
    return k, v


def label(point):
    """
	Predict the label for a given datapoint using MLlib Kmeans 'clusters' object
	"""
    return clusters.predict(point)


if __name__ == "__main__":
    # Load from file
    spikes = load('test_data.npy')

    # Reformat dataset
    spikes[:, 0] -= 1  # To establish zero-indexing
    spikes[:, 1] *= 1e3  # Convert spike times into milliseconds
    train_period = 1000.  # Period of spike train, for plotting later
    n_trains = 105  # Number of spike trains
    spikes[:, [0, 1]] = spikes[:, [1, 0]]  # Switch column order

    # Create binary vectors for each neuron. This is the Spark friendly format
    trains = zeros((n_trains, train_period))
    spk = spikes.astype(int)
    trains[spk[:, 1], spk[:, 0]] = 1

    # Keep only those trains that have spikes (all here, though how to keep track of un-clustered trains?)
    active = trains.any(axis=1)
    act_trains = trains[active, :]

    # Convert the data into a spark RDD
    # sc = PySparkShell.start(appName="community")
    sc = SparkContext("local", "community_detection")
    keys = arange(1, shape(trains)[0] + 1)
    spks = sc.parallelize(zip(keys, trains))

    # Compute ISIs of all trains in parallel
    ISIs = spks.mapValues(compISIs)

    # Collect ISIs and compute max and min
    ISIstats = ISIs.map(lambda (k, v): (k, [median(v), min(v)] )).values().collect()

    x, y = zip(*ISIstats)
    bin_max = median(x)
    bin_min = min(y)
    # print bin_max
    # print bin_min

    # Compute timescale distribution - check how Mark does it
    # Set up binwidths
    t_reps = 7
    bins = linspace(bin_min, bin_max, t_reps)

    # For binless version, divide binwidth by sqrt(12), and remove timescales of < 1ms
    timescales = bins / sqrt(12)
    timescales = timescales[timescales >= 1.]
    print timescales[0]

    # Convolution. This is the start of the main loop.
    # i = 1
    spks_conv = spks.map(lambda (k,v): (k, convolve_series(v, timescales[0])))

    # Collect the convoluted spike trains, then ship the matrix to each node with the compute_similarity function.
    # This method may be cheaper (depending on matrix size i.e. communication costs) this way, than for each node to have to query one another all the time.
    # NOTE: Could also try broadcasting SC without collect(). How is this done? Is this handled more efficiently?
    SC = spks_conv.collect()
    SC = zip(*SC)[1]
    SC = array(SC)
    # plt.plot(SC)
    # plt.show()

    # Compute similarity matrix
    A = spks_conv.map(lambda (k,v): (k, compute_similarity(v, SC)))
    # Set diagonal elements to zero
    A = A.map(lambda (k,v): zero_diag(k,v))

    # Convert A into a Thunder RowMatrix for more flexible manipulation
    A_rm = matrices.RowMatrix(A)

    # Compute Expected matrix
    # Compute sum(A_rm, axis=0) and sum(A_rm), to get the outer product of the two axis sums without having to sum across pipelined RDD rows
    r = A_rm.rows().sum()
    n = sum(r)
    outer_A = A_rm.rdd.mapValues(lambda x: sum(x)*r)
    E = outer_A.mapValues(lambda x: x/n)

    # Finally, subtract E from A to generate B, the modularity matrix
    E_rm = matrices.RowMatrix(E)
    B = A_rm.elementwise(E_rm,'minus')

    # Eigendecomposition with Thunder's SVD
    # need a way of keeping as many eigenvalues as needed
    svd = SVD(k=20, method="direct")
    svd.calc(B)
    n_units, n_features = svd.v.T.shape
    data = svd.u

    # Collect a local copy for computing Q. Should be done with a reduce step instead
    b_local = B.collect()

    # Cluster with MLlib Kmeans
    Q = 0
    L = empty((n_features, n_units))

    # Build the model (cluster the data)
    for i in range(n_features):
        K = i+2
        clusters = KMeans.train(data.values(), K, maxIterations=100, runs=100, initializationMode="k-means||")

        Labels = data.mapValues(lambda point: label(point)).collect()

        # Construct membership matrix
        S = zeros((n_units, K))
        for n in xrange(n_units):
            S[n, Labels[n][1]] = 1

        # Compute clustering quality (eq. 11 in the paper)
        q = matrix.trace(dot(dot(S.T, b_local), S))

        # Save best result
        if q > Q:
            Q, L, C = q, Labels, clusters.centers

    print Q

    # Assign group labels to original spike train and plot
    n_clusters = array(C).shape[0]
    outsp, outn = [[] for x in xrange(n_clusters)], zeros(n_clusters)
    for i in xrange(array(L).shape[0]):
        label = L[i][1]
        for j in xrange(spikes.shape[0]):
            if spikes[j,1] == i:
                spt = spikes[j,0]
                spidx = outn[label]
                outsp[label].append((spt, spidx))
        outn[label] += 1


    # Plot results
    fig, ax = plt.subplots(2,1,figsize=(12, 8))
    colors = cycle(['b', 'g', 'r', 'c', '', 'y', 'k', 'w'])

    ax[0].plot(spikes[:,0],spikes[:,1],'.',ms=4.5)
    ax[0].set_xlabel('Time (ms)')
    ax[0].set_ylabel('Cell ID')



    for i in xrange(n_clusters):
        outsp[i] = array(outsp[i])
        ax[1].scatter(outsp[i][:,0], outsp[i][:,1] + sum(outn[:i]),
            color=next(colors), s=2.5)
        ax[1].set_xlabel('Time (ms)')
        ax[1].set_ylabel('Cell ID')

    plt.xlim([0,1000])
    plt.ylim([0,120])
    plt.show()



