%%% script to run binned and binnless clustering methods, and plot results
%%% This script demonstrates finding spike-train groups based on
%%% millisecond-scale correlations in spike-times. 
% Mark Humphries 02/06/10
close all; clear all

load test_data
% loads a single data-set from Fellous et al (2004) test data. This set has
% 3 different spike-train responses, each response lasting 1 second. Each
% response is 4-6 spikes at fixed times. Each response occurs for 35 trials (total of 105 trials). 
% On any trial, there is independent added noise from 3 sources:
% (1) 15% chance of each fixed spike missing
% (2) 2 extra spikes added at random times
% (3) A Gaussian-derived jitter of each remaing fixed spike, with s.d. of 1ms
% Note: this is the least "noisy" test set, used here for clarity in the
% end plots

% rand('state',1); rand('state',1);
% data parameters
ntrials = 105;
T = [0,1];  % all trials are 1 second

% binned options
binopts.Dmeth = 'Hamming';    % Hamming similarity
binopts.modopts = {'',20};  % use Euclidean distance for k-means (default), repeat 20 times

% binless options
step = 0.001;   % binless window step in seconds
binlessopts.Dmeth = 'cosine';   % use cosine of angle as similarity
binlessopts.modopts = {'',20};  % use Euclidean distance for k-means (default), repeat 20 times

%% run analysis
% get spike-train stats
ISIs = []; Hz = zeros(ntrials,1);
for i = 1:ntrials
    currix = find(spks(:,1) == i);
    ts = spks(currix,2);
    ISIs = [ISIs; diff(ts)];
    Hz(i) = 1 ./ mean(diff(ts));
end
mISI = mean(ISIs); medISI = median(ISIs);
[eisis,x] = ecdf(ISIs); % empricial CDF of ISIs...
isi99 = x(find(eisis <= 0.01,1,'last'));    % use 1% ISI dist as smallest bin

%% find clusters - using bins
% dynamic bin-size
mInterval = min(mISI,medISI);           % centre interval
binsize = linspace(isi99,mInterval,7);  % try 7 different binsizes
[G,GS,Ngrps,Qm] = cluster_spike_data_bins(spks,[1:ntrials]',binsize,T,binopts);

ix = find(Ngrps == 3); ix = ix(1);  % choose binsize with (possibly) correct grouping
hc = plot_clusters(spks,G{ix},Ngrps(ix),T,'123');   % look at found groups

%% find clusters - using binless
binlessopts.BLpars = binsize / sqrt(12);   % following Kruskal et al (2007)
[Gbinless,GSbinless,Ngrps_binless,Qm_binless] = cluster_spike_data_binless(spks,[1:ntrials]',T,step,binlessopts);

ix = find(Ngrps_binless==3); ix = ix(1);
hc = plot_clusters(spks,Gbinless{ix},Ngrps_binless(ix),T,'123');   % look at found groups


%% find clusters - using bins and permutation test: constructs surrogate
% spike-train data to find significantly similar pairs, and constructs
% network using only those links [same options are required for running
% this with "binless" function]

binopts.blnGrph = 1; 
binopts.alpha = 0.05; 
binopts.npermutes = 100;

[Gp,GSp,Ngrpsp,Qmp] = cluster_spike_data_bins(spks,[1:ntrials]',binsize,T,binopts);

pix = find(Ngrps == 3); pix = pix(1);  % choose binsize with (possibly) correct grouping
hc = plot_clusters(spks,G{pix},Ngrps(pix),T,'123');   % look at found groups

%% using binless results, further examples of plotting
% n.b. here we plot only the detected groups

% (1) specifying colours
cmap = [1 0 0; 0 0.5 0; 0 0 1]; % plot 3 groups as red, dark green, and blue   
hc2 = plot_clusters(spks,Gbinless{ix},Ngrps_binless(ix),T,'3',[],cmap);
title('Example use of custom color map')

% (2) specifying order of plotting (for grouped plot), by firing rate; 
%     plots groups such that top spike-train is most active, bottom spike-train is least active 
%     n.b. vector Hz containing firing rates computed above
hc2 = plot_clusters(spks,Gbinless{ix},Ngrps_binless(ix),T,'3',Hz,cmap);
title('Example use of intra-group plotting by order: firing rate')

%% Advanced plotting: specifiying order of plotting by intra-group similarity
% n.b. these plots are of limited value in this example where the
% spike-trains are all highly structured; however, these forms of plots are
% very useful when examining multi-neuron data

% (1) order spike-trains within groups by their total similarity to other
% members of the group (top spike-train is most similar, bottom spike-train
% is least similar)

% use new blnS option to run function without clustering, and return
% similarity matrix only for specified Gaussian width
allwdths = binlessopts.BLpars;
binlessopts.BLpars = binlessopts.BLpars(ix); 
binlessopts.blnS = 1; [Gx,GSx,Ngrpsx,Qx,Sx] = cluster_spike_data_binless(spks,[1:ntrials]',T,step,binlessopts);

tS = sum(Sx{1});  % total similarity with rest of network
Sin = zeros(ntrials,1); Sout = zeros(ntrials,1); Sgrp = zeros(Ngrps_binless(ix),1); 
thisG = Gbinless{ix};
for j = 1:ntrials
    % compute each spike-train's similarity within their group, and outside
    % of their group...
    thisgrp = thisG(j,2);
    ingrp = find(thisG(:,2) == thisgrp); 
    outgrp = find(thisG(:,2) ~= thisgrp);     
    Sin(j) = sum(Sx{1}(j,ingrp));   % total similarity of current spike-train with members of its group
    Sout(j) = sum(Sx{1}(j,outgrp)); % total similarity of current spike-train with all spike-trains not in its groups
end
hc2 = plot_clusters(spks,Gbinless{ix},Ngrps_binless(ix),T,'3',Sin,cmap);
title('Example use of intra-group plotting by order: total similarity')

% (2) same again, but now we also re-map group membership indices so that groups
% are also plotted by mean total intra-group similarity (so that most
% similar group at the top, least similar group at the bottom)

for j = 1:Ngrps_binless(ix)
    ingrp = find(thisG(:,2) == j); % find all group members 
    Sgrp(j) = sum(Sin(ingrp))/numel(ingrp); % mean total intra-group similarity
end
[srted,s_idx] = sort(Sgrp);    % low-to-high
tempG = thisG; 
for j = 1:Ngrps_binless(ix)  thisG(tempG(:,2) == s_idx(j),2) = j; end  % remap group membership indices

hc3 = plot_clusters(spks,thisG,Ngrps_binless(ix),T,'3',Sin,cmap);
title('Plotting groups by mean total similarity, plotting within groups by total similarity')

