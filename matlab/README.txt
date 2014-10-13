%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spike-Train Communities Toolbox v1.2				   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Put this folder of code somewhere, and add to your MATLAB path.

All functions have extensive help documentation available from the MATLAB prompt.

USAGE:
To run the whole process (construct similarity matrix and run algorithm for every specified time-scale) simply pass time-series data to the top-level functions. 

If similarity matrices are already available (e.g. have already been computed from spike-train data using alternative correlation metrics, or from imaging time-series) then use allevsplit.m to directly run the community detection algorithm

For interpretation of output, and further controls using shuffled data, see Humphries (2011) J Neurosci

FILE LIST:

example_script.m: shows how to call the top-level functions with time-series data, and how to visualise the results in various ways

Top-level functions :
cluster_spike_data_bins.m: takes time-series in array format (event ID, time-stamp), and specification of similarity measure and time-scale(s) for bins over which these are to be computed. It then computes the similarity matrix from the time-series at each specified time-scale, and runs the community-detection algorithm on each of those matrices. Returns the group structure and modularity scores found at each time-scale. Has extensive options for setting up the matrices, including setting a permutation-test derived threshold for defining a binary adjacency matrix, should this be required.

cluster_spike_data_binless.m : as above, but convolves the time-series with the specified kernel (Gaussian, exponential) at the specified range of time-scales, and computes the specified similarity metric.

Visualisation:
plot_clusters.m : takes the time-series array and group structure, and plots the rasters in group order, colour-coded.	

Algorithm functions:
allevsplit.m : takes the matrix defining the network (whether adjacency matrix, or similarity matrix), and runs the eigenspectra-based clustering algorithm from Humphries (2011) J Neurosci.
Returns the grouping and modularity score.

multileadevsplit.m : implements the leading-eigenvector methods from Newman (2006) Phys Rev E, which iteratively bifurcates the network into sub-networks. Works well on standard undirected networks (adjacency matrices). 

expectedA.m : constructs the matrix (P) for the expected network, called by both allevsplit.m and multileadevsplit.m

refinesplit.m : [optional step for multileadevsplit.m] at each step of the algorithm swaps each node into the other group, and recomputes modularity. Retains all swaps that increase modularity. 


Helper functions:
discreteinvrnd.m : generate random numbers drawn from an arbitrary discrete probabiltiy distribution

MIpartitions : compute normalised mutual information between two groupings of the same network. Useful for comparing group structure of the same neurons over different time-points (see e.g. Fig 7 of Humphries(2011) J Neurosci). 

shuffle_intervals.m: takes spike-train data-set, and shuffles individual spike-train time-series tom construct control data-set. Control data-sets can be passed to the top-level function to generate null model distribution for modularity

spike_train_from_times.m : construct binary vector of time-series from compressed representation in array format (eventID,time).  