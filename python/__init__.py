"""
Implementation of Spike-Train Communities algorithm (Humphries, 2011).

`find_assemblies()' is the `master' function of this file.

Thomas Sharp, 2013
thomas.sharp@riken.jp
"""
import numpy
import scipy.cluster
import scipy.spatial
import scipy.stats

DEBUG = False



def cluster_graph(A, P, reps, n_clusters=None):
    """
    Performs clustering on the measure of pairwise similarity between convolved
    spike trains. See the paper for details.
    
    :param numpy.ndarray A:
        pairwise-similarity measure between convolved spike trains.
    :param numpy.ndarray P:
        null-model.
    :param int reps:
        number of times to repeat the clustering.
    
    :param int n_clusters:
        (optional) fixed number of clusters.
    """
    # Compute modularity matrix
    B = A - P
    # Take eigenvectors of positive eigenvalues as clustering features
    eigenvalues, eigenvectors = numpy.linalg.eigh(B)
    features = eigenvectors[:,eigenvalues>0] # One feature per column
    n_units, n_features = features.shape
    # Run clustering
    Q, L, centres = 0, None, None
    if not n_clusters: n_clusters = n_features + 1
    if n_clusters > n_features + 1:
        raise Exception('There cannot be more clusters than features.')
    L = numpy.empty((n_features, n_units))
    for clusters in range(2, n_clusters + 1):
        if DEBUG: print "    ... on %d clusters" % clusters
        # Repeat 'reps' times for stability
        for j in range(reps):
            if DEBUG: print "        ... k-means run %d of %d" % (j+1, reps)
            # Generate initial cluster centres
            init = initialise_centres(features, clusters)
            # Prepare values for clustering and run kmeans
            whitened = scipy.cluster.vq.whiten(features)
            test_centres, labels = scipy.cluster.vq.kmeans2(whitened, init)
            # Construct membership matrix
            S = numpy.zeros((n_units, clusters))
            for k in range(n_units):
                S[k, labels[k]] = 1
            # Compute clustering quality (eq. 11 in the paper)
            q = numpy.matrix.trace(numpy.dot(numpy.dot(S.T, B), S))
            # Save result
            if q > Q:
                Q, L, centres = q, labels, test_centres
    return Q, L, centres
    
    
def compute_similarity(trains):
    """
    Computes pairwise similarity between all spike trains and produces a
    null-model network for comparison.

    :param numpy.ndarray trains:
        2d array with neurons on columns, time on rows and spikes represented
        by ones.
    
    :returns:
        tuple of two numpy.ndarray instances containing the similarity measure
        and a null-model for comparison.
    """
    # Find similarity between continuous spike-train vectors
    A = scipy.spatial.distance.pdist(trains, 'cosine')
    A = scipy.spatial.distance.squareform(1 - A) # similarity = 1 - distance
    # Make null-model 'expected' network (eq. 4 in the paper)
    P = numpy.outer(numpy.sum(A, axis=0), numpy.sum(A, axis=1)) / numpy.sum(A)

    return A, P


def compute_timescales(spikes, n, reps):
    """
    Computes candidate timescales for spike train comparison based on interspike
    intervals (ISIs).

    :param numpy.ndarray spikes:
        2d array with neurons on columns, time on rows and spikes represented
        by ones.
    :param int n:
        number of spike trains in the input.
    :param int reps:
        number of timescales to return.

    :returns:
        numpy.ndarray containing candidate timescales for the analysis.
    """
    # Gather all ISIs
    ISIs = list()
    for i in range(n):
        idx = spikes[:,1] == i
        spk = spikes[idx,0]
        isi = numpy.diff(spk)
        ISIs.extend(isi)
    ISIs = numpy.sort(ISIs)
    # Compute empirical CDF of ISIs; take 1%-point as smallest bin
    ecdf = numpy.arange(ISIs.size) / float(ISIs.size)
    idx = numpy.abs(ecdf - 0.01).argmin()
    bin_min = ISIs[idx]
    # Produce candidate timescales, according to ISI statistics
    mean, median = numpy.mean(ISIs), numpy.median(ISIs)
    bin_max = numpy.min([mean, median])
    bins = numpy.linspace(bin_min, bin_max, reps)
    # Divide according to source paper, and remove sigmas of < 1ms
    sigmas = bins / numpy.sqrt(12)
    sigmas = sigmas[sigmas >= 1.]

    return sigmas


def convolve_trains(trains, sigma):
    """
    Convolves spike trains with a Gaussian to produce pseudocontinuous signals.

    :param numpy.ndarray trains:
        2d array with neurons on columns, time on rows and spikes represented
        by ones.
    :param int sigma:
        timescale (milliseconds) on which to build filter.

    :returns:
        numpy.ndarray of the same shape as the trains input, but containing
        continuous values representing the convolved trains.
    """
    # Find train parameters
    n, t = trains.shape
    # Generate the filter and output arrays
    sigma = numpy.int(numpy.round(sigma))
    window = numpy.arange(-5*sigma, +5*sigma)
    gaussian = scipy.stats.norm.pdf(window, scale=sigma)
    conv_trains = numpy.empty((n, t+10*sigma-1))
    # Make a continuous signal for each train by convolution
    for i in range(n):
        conv_trains[i,:] = numpy.convolve(trains[i,:], gaussian)

    return conv_trains


def initialise_centres(features, n_clusters):
    """
    Finds initial cluster centres for K-means clustering according to the
    K-means++ algorithm (Arthur and Vassilvitskii, 2007).

    :param numpy.ndarray features:
        data on which the clustering algorithm will be run, with one sample
        on each row and one feature per column.
    :param int n_clusters:
        number of clusters to initialise for.

    :returns:
        numpy.ndarray of cluster centres, with shape (n, f) for n clusters
        and f feature dimensions.
    """
    # Extract meta-information and set up return array
    n_units, n_features = features.shape
    centres = numpy.empty((n_clusters, n_features))
    # Choose an initial centre at random
    init = numpy.random.randint(0, n_units)
    centres[0,:] = features[init,:]
    # Choose successive centres with prob. func. of distance
    for i in range(1, n_clusters):
        d = numpy.empty(n_units)
        # Find distance from each point to nearest centre
        for j in range(n_units):
            p_tmp = features[j,:] # Current point
            c_tmp = centres[:i,:] # Existing centres
            d_tmp = numpy.sqrt(numpy.sum((c_tmp - p_tmp) ** 2, axis=1))
            d[j] = numpy.min(d_tmp) # Distance to nearest centre
        # Select a new center with probability proportional to d^2
        d_prob = numpy.cumsum(d**2 / numpy.sum(d**2))
        rnd = numpy.random.uniform(0,1)
        idx = d_prob.searchsorted(rnd)
        centres[i,:] = features[idx,:]

    return centres


def find_assemblies(n_trains, train_period, spikes, t_reps=7, k_reps=10,
    n_clusters=None):
    """
    Master function of this module. Takes spike trains and returns clustering
    quality, the cluster label of each spike train, the cluster centers and
    the timescale upon which the best clustering was found.
    
    :param int n_trains:
        number of spike trains in the input.
    
    :param int train_period:
        length of the spike trains.
    
    :param numpy.ndarray spikes:
        spikes in two columns corresponding to firing times (in milliseconds)
        and cell IDs.
    
    :param int t_reps:
        number of timescales on which to run the algorithm.
    
    :param int k_reps:
        number of times to iterate k-means.
    
    :param int n_clusters:
        (optional) fixed number of clusters.
    
    :returns:
        A 4-tuple containing:
        
            Q (int): quality of the clustering (see paper eq. 11),
            
            L_out (numpy.ndarray): cluster labels for each of the n_trains,
            
            C (numpy.ndarray): cluster centres, with shape (n, f) for n clusters
                and f feature dimensions.
            
            T (int): timescale of the best clustering (milliseconds)
    """
    Q, L, C, T = 0, None, None, None
    trains = format_spikes(spikes, n_trains, train_period)
    # Filter trains with no spikes in them
    active = trains.any(axis=1)
    act_trains = trains[active,:]
    # Run the clustering algorithm for each candidate timescale
    timescales = compute_timescales(spikes, n_trains, t_reps)
    for i in range(timescales.size):
        if DEBUG: print "Running at timescale %dms" % int(timescales[i])
        # See the paper for detailed algorithm details
        conv_trains = convolve_trains(act_trains, timescales[i])
        A, P = compute_similarity(conv_trains)
        q, l, c = cluster_graph(A, P, k_reps, n_clusters)
        if q > Q:
            Q, L, C, T = q, l, c, timescales[i]
    # Assign labels to the "no-spike" trains that weren't filtered out
    L_out = numpy.ones(n_trains, dtype=int) * -1
    L_out[active] = L

    return Q, L_out, C, T
    

def format_spikes(spikes, n, t):
    """
    Convert spike array from sparse to dense encoding.

    :param numpy.ndarray spikes:
        spikes in two columns corresponding to firing times in milliseconds and
        cell IDs in unique integers.
    :param int n:
        number of spike trains in the input.
    :param int t:
        number of timesteps in the input.

    :returns:
        2d array with neurons on columns, time on rows and spikes represented
        by ones.
    """
    # Set up output and sanitise spike type
    trains = numpy.zeros((n, t))
    spk = spikes.astype(numpy.int)
    # Write spikes to binary matrix
    trains[spk[:,1], spk[:,0]] = 1

    return trains
