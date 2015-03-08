"""
Functions used by community_detection.py 
Written by Mathew H. Evans 2015. Based on 'Spike train community detection'. Humphries, 2011.
"""
from thunder.utils import pack, getdims, matrices
from thunder.factorization import PCA, SVD
# from thunder.clustering.kmeans import KMeans, KMeansModel
from thunder.viz import Colorize
import numpy as np
import itertools
import pylab
import scipy.spatial
import scipy.stats
import scipy.cluster
import time 
from pyspark.mllib.clustering import KMeans
from math import sqrt




def compISIs(self):
    """
    Compute ISIs for each spike train
    """
    x = np.where(np.array(self)>=1)
    
    return np.diff(x)


def convolve(self, timescale):
	"""
	Convolve each spike train with a Gaussian of width 'timescale'. 
	Returns the convolved spike train
	"""
	import scipy.stats

	# Generate the filter and output arrays
	sigma = np.int(np.ceil(timescale))
	window = np.arange(-5*sigma, +5*sigma)
	gaussian = scipy.stats.norm.pdf(window, scale=sigma)
	conv_trains = np.empty(0) #(0,(np.array(self).shape[1])+10*sigma-1))

	# Make a continuous signal for each train by convolution

	conv_train = np.convolve(np.array(self), gaussian)

	return conv_train

def compute_similarity(self,SC):
    """
    Computes pairwise similarity between all spike trains using scipy cosine similarity
    """
    import scipy.spatial
    
    x = self #np.transpose(self)
    
    # Find similarity between continuous spike-train vectors
    A = np.empty([SC.shape[0]])
    
    for i in range(SC.shape[0]):
        A[i] = scipy.spatial.distance.pdist([x,SC[i]], 'cosine')
        
    # Similarity = 1 - distance
    A = 1 - A

    return A

def zero_diag(k,v): 
	"""
	Remove diagonal elements of a matrix 
	"""
	v[k-1] = 0
	return k,v

def label(point):
	"""
	Predict the label for a given datapoint using MLlib Kmeans 'clusters' object
	"""
	return clusters.predict(point)