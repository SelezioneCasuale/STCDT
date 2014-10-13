function [G,grpsizes,ngrps,Beta,varargout] = cluster_spike_data_binless(spkdata,Didxs,T,bin,varargin)

% CLUSTER_SPIKE_DATA_BINLESS cluster spike-data using binless-based metrics
% [G,GI,GS,C] = CLUSTER_SPIKE_DATA_BINLESS(S,I,T,B,OPTS) clusters spike trains into
% groups, based on similarities between pairs of "binless" spike train representations.
% Similarities between all pairs are stored in similarity matrix, which forms
% the basis for the clustering.
%
% S is a 2-column vector: first column is ID stamps; second column is
% time-stamps. IDs can be either for events (e.g. stimulus presentations
% for single unit data; time-stamps are then relative to stimulus presentation)
% or neurons (for multi-unit data; time-stamps are then relative to recording period).
%
% I is a column vector of all possible ID stamps - note that they need not all appear in S
% because the ID event may have caused no spikes or been a silent neuron!
% If a neuron is silent, it will not appear in the groupings: its ID stamp
% will be omitted from G.
%
% T is 2-element vector of [start end] for all time-stamps (e.g. if using zero as
% time-stamp for stimulus presentation: for 2 seconds of data starting at stimulus presentation T = [0 2]; for 2
% seconds starting 0.5 second before stimulus T = [-0.5 1.5];
%
% B is the spike-time quantisation value - default is 1ms (B=0.001)...
%
% ... = CLUSTER_SPIKE_DATA_BINS(...,OPTS) sets all options contained in the OPTS structure - omitting a field will use the default option, indicated
% by the parentheses:
%       OPTS.BLmeth = {('Gaussian') | 'exponential'}
%               Convolve each spike with a Gaussian (Fellous et al, 2004; used in Humphries, 2010); or
%               forward exponential function (van Rossum, 2001; included for future development of algorithm)
%       OPTS.BLpars : array of values for the binless metric parameter - each value will be tested,
%               'Gaussian' - standard deviation of Gaussian in seconds (default: 0.01)
%               'exponential' - decay constant of exponential in seconds
%       OPTS.Dmeth = {{'cosine'}| 'corrcoef' | 'corr'}
%           sets the comparison-between-pairs method to: cosine of angle between vectors (Fellous et al 2004; used in Humphries, 2010);
%           rectified correlation coefficient (corrcoef or corr)
%       OPTS.blnGraph = {(0)|1} : converts similarity matrix into adjacency
%           matrix using permutation test to find p-value for similarity
%           between each pair of spike-trains
%       OPTS.alpha (0.05) : sets the significance threshold for graph conversion (has no effect otherwise)
%           pairs of spike-trains with p < alpha are considered significantly similar,
%           so are assigned a 1 in the adjacency matrix, everything else is assigned a 0.
%       OPTS.npermutes (100) : number of permutations (shuffles of each
%           train) to get p-values; if sufficient computational power is
%           available, then 1000 is a good starting point
%       OPTS.nlimit (6) : the minimum number of retained entries from graph conversion required to
%           proceed with the clustering analysis (default is 6: two groups
%           of 3 nodes each) [has no effect if blnGraph = 0].
%       OPTS.blnAll = {0|(1)} - uses the all positive eigenvector method to cluster the spike trains
%           (function ALLEVSPLIT) as in Humphries, 2010: this works best on full correlation
%           matrices (i.e. blnGraph=0); setting blnAll = 0 uses the leading
%           eigenvector method instead (function MULTILEADEVSPLIT), which
%           works better on adjacency matrices.
%       OPTS.modopts ({}) - passes anything in this cell array to the corresponding optional
%           arguments of the clustering functions (ALLEVSPLIT or MULTILEADEVSPLIT) - see
%           their help for details. For example, for ALLEVSPLIT, the first cell will be
%           passed into the first optional argument (a string),
%           the second cell into the second optional argument (a number)
%       OPTS.blnS = {(0)|1} - if set to 1, does not do the clustering algorithm, but only computes
%           the similarity matrices for each value in opts.BLpar: useful for just creating a set of such matrices for further processing.
%           [Default is 0]
%
% Returns: cell array of all group structure G, one cell per binsize
% tested; each with a two-column vector: first column is ID stamp, second column
% is group membership (integer value) - note that there is no guarantee all IDs will be assigned a group.
% A cell array of group sizes GS; and measure of grouping quality C: for the ALLEVSPLIT algorithm (Humphries, 2010),
% this just returns maximum modularity Q; for MULTIEVSPLIT, this returns metric in (Humphries et al 2009).
%
% [...,Sxy,BD,Axy] = CLUSTER_SPIKE_DATA_BINLESS(...) are optional outputs useful for further
% processing: Sxy is a cell array of the similarity matrices; BD is a cell array of
% matrices containing the convolved spike-train vectors - each column is one
% spike train; Axy is a cell array of the adjacency matrices resulting from setting blnGraph=1 (is empty otherwise).
% For each, there is one cell per tested binless metric parameter.
%
% Notes:
% (1) This is version 2 of this function. The main change is the option to construct an unweighted, undirected graph
%     from the data, using a permutation test to define significant
%     similarity, and thus to define which spike-train pairs are linked in
%     the graph. Can only do this with relatively small (N of a few
%     hundred) data-sets, otherwise will rapidly run out of memory and
%     life-span
%
% (2) The novel algorithm described in Humphries (2010) corresponds to the
% default choices in this function. Other options are included for further
% exploration of applications of this analysis technique.
%
% (3) Choice of Gaussian widths: typically, we choose a set of discrete binsizes based on
%  ISI statistics of spike-train data-set, then follow Kruskal et al (2007)
%  and convert to Gaussian widths by s = w / sqrt(12)
%
% (4) If an ID stamp has no associated spikes (e.g. there was no response
% to the stimulation or an error in recording) then the train with that ID
% is omitted from the analysis: the groupings in G will not contain that ID
% stamp.
%
% References:
% (1) Humphries, M. D. (2011) Spike-train communities: finding groups of similar
% spike trains J Neurosci 31, 2321-2336
% (2) Humphries, Wood & Gurney (2009) Dopamine-modulated dynamic cell assemblies generated by the GABAergic striatal
% microcircuit, Neural Networks, 22, 1174-1188.
% (3) van Rossum M. C. (2001) A novel spike distance. Neural Comput, 13, 751-763
% (4) Fellous, J. M.; Tiesinga, P. H.; Thomas, P. J. & Sejnowski, T. J.
%       (2004) Discovering spike patterns in neuronal responses J Neurosci, 24, 2989-3001
% (5) Kruskal, P. B., Stanis, J. J., McNaughton, B. L., Thomas, P. J. (2007). A binless correlation measure reduces the variability
% of memory reactivation estimates. Stat Med 26: 3997–4008.
%
% Mark Humphries 12/6/2011

[r c] = size(Didxs);
if r==1 Didxs = Didxs'; end   % make column vector

% defaults
bin = 0.001;    % 1ms bins
opts.BLmeth = 'Gaussian';   % use Gaussian window around spikes
opts.BLpars = 0.01;        % std dev is 10 ms
opts.Dmeth = 'cosine';     % use correlation coefficient
opts.blnGraph = 0;   % don't make graph
opts.nlimit = 6;     % minimum number of retained nodes for graph...
opts.alpha = 0.05;    % signficance threshold for retaining in graph
opts.npermutes = 100;   % default number of permutations

% opts.tau = [0.01:0.01:0.3];      % tested interval of slopes of dynamic range amplification function (see Fellous et al 2004)

opts.blnAll = 1;     % use all eigenvector method
opts.modopts = '';   % set no options for the all eigenvector method...
opts.blnS = 0;       % run clustering algorithm by default

if nargin >= 5
    if isstruct(opts)
        tempopts = varargin{1};
        fnames = fieldnames(tempopts);
        for i = 1:length(fnames)
            opts = setfield(opts,fnames{i},getfield(tempopts,fnames{i}));
        end
    end
end

% storage
Rn = cell(numel(opts.BLpars),1);        % retained neurons...
grpsizes = cell(numel(opts.BLpars),1);
ngrps = zeros(numel(opts.BLpars,1));
spkfcn = cell(numel(opts.BLpars),1);

bins = T(1)+bin/2:bin:T(2);

%% analyse data
for loop = 1:numel(opts.BLpars)
    
    % set up convolution window if using...
    sig = opts.BLpars(loop)/bin;
    switch opts.BLmeth
        case 'Gaussian'
            x = [-5*sig:1:5*sig]';  % x-axis values of the discretely sampled Gaussian, out to 5xSD
            h = (1/(sqrt(2*pi*sig^2)))*exp(-((x.^2*(1/(2*sig^2))))); % y-axis values of the Gaussian
            shiftbase = floor(length(x)/2); % was floor....
            % keyboard
        case 'exponential';
            x = [0:1:10*sig];   % spread to 10 times the time constant
            h = exp(-x/sig);
            shiftbase = length(x);
    end
    
    % convolve window with data-spike trains
    [spkfcn{loop},idxs] = convolve_spiketrains(spkdata,h,shiftbase,Didxs,bins,bin,T,opts);
    Nidxs = numel(idxs);
    
    %% now compute selected distance between those functions
    [Sxy{loop}] = constructS(spkfcn{loop},Nidxs,opts);
    
    %% if constructing adjacency matrix, then run permutation tests and
    %% find graph structure...
    if opts.blnGraph
        % run permutation test
        for ctrl = 1:opts.npermutes
            % get shuffled trains
            ctrlspkts = shuffle_intervals(spkdata,T,bin);
            % convolve window with control spike trains
            [ctrlspkfcn,ctrlidxs] = convolve_spiketrains(ctrlspkts,h,shiftbase,Didxs,bins,bin,T,opts);
            
            %% now compute selected distance between those functions
            ctrlSxy{ctrl} = constructS(ctrlspkfcn,Nidxs,opts);
        end
        
        % for each entry in data S, find significant similarities, and
        % construct adjacency matrix
        Ac = zeros(Nidxs);
        for cc = 2:Nidxs
            for rr = 1:cc-1
                thisSdist = zeros(opts.npermutes,1);
                for btch = 1:opts.npermutes
                    thisSdist(btch) = ctrlSxy{btch}(rr,cc);
                end
                thisSdist = sort(thisSdist);
                Ac(rr,cc) = Sxy{loop}(rr,cc) > prctile(thisSdist,100*(1-opts.alpha));   % find if data similarity sig greater than shuffled
                Ac(cc,rr) = Ac(rr,cc);
            end
        end
        
        % strip out un-connected nodes
        dgr = sum(Ac);
        delN = find(dgr < 2);
        Rn{loop} = setdiff(1:Nidxs,delN);    %% index back into Rn to recover original ID.
        Acs{loop} = Ac(Rn{loop},Rn{loop});
        
    else
        Acs{loop} = Sxy{loop};
        Rn{loop} = 1:Nidxs;
    end
    
    % keyboard
    % do modularity on that graph
    edges(loop) = sum(sum(Acs{loop}~=0)); nodes(loop) = numel(Rn{loop});
    if ~opts.blnS & nodes(loop) >= opts.nlimit & edges(loop) > log(nodes(loop))
        % only group if: (a) asked to do so;
        % (b) there are enough nodes and
        % (c) the number of edges ensures a likely fully-connected graph
        
        if opts.blnAll
            if numel(opts.modopts) == 1
                [Sc,Vpos,B,Q] = allevsplit(Acs{loop},opts.modopts{1});
            elseif numel(opts.modopts) == 2
                [Sc,Vpos,B,Q] = allevsplit(Acs{loop},opts.modopts{1},opts.modopts{2});
            else
                [Sc,Vpos,B,Q] = allevsplit(Acs{loop});
            end
        else
            if numel(opts.modopts) > 1
                error('Too many options passed to MULTILEADEVSPLIT')
            else
                [Sc,Uc] = multileadevsplit(Acs{loop},opts.modopts);
            end
        end
        
        ngrps(loop) = max(Sc);
        % get sizes
        siz = [];
        for i = 1:ngrps(loop)
            siz = [siz numel(find(Sc == i))];
        end
        grpsizes{loop} = siz;
    else
        Sc = zeros(numel(Rn{loop}),1);  % nothing to group!
        ngrps(loop) = 0;
        grpsizes{loop} = [];
        Q = 0;
        B = zeros(Nidxs);
    end
    
    % remap from retained index count to ID stamps...
    G{loop} = [idxs(Rn{loop}) Sc];  % [ID stamp; Group membership]
    
    if opts.blnAll
        % then return max Q as the measure of grouping
        maxQ =  max(max(Q)); if isempty(maxQ) maxQ = 0; end
        Beta(loop) = maxQ;
        
    else
        %otherwise if used MULTILEADEVSPLIT, must use some other metric:
        %this uses the NN (2009) paper metric
        nzSxy = Sxy{loop}(Sxy{loop} > 0);
        Beta(loop) = sum(grpsizes{loop})/Nidxs * ngrps(loop) * abs(min(nzSxy) - median(nzSxy));
    end
end


varargout{1} = Sxy;
varargout{2} = spkfcn;
varargout{3} = Acs;

function [spkfcn,idxs] = convolve_spiketrains(spkdata,h,shiftbase,Didxs,bins,bin,T,opts)
Nidxs = numel(Didxs);

%% go round and compute spike-train binless functions
spkfcn = zeros(numel(bins),Nidxs);

for j = 1:Nidxs
    currix = find(spkdata(:,1) == Didxs(j));
    nspikes(j) = numel(currix);
    [spk,bts] = spike_train_from_times(spkdata(currix,2),bin,T);
    if nspikes(j) > 0   % only bother doing convolution if there's something to convolve!!
        switch opts.BLmeth
            case 'Gaussian'
                try
                    y = conv(h,spk);
                    [r c] = size(y); if c>1 y = y'; end  % for reasons best known to Matlab, certain convolutions will return this as a row vector rather than a column vector
                    shifty = y(shiftbase+1:end-shiftbase);   % shift convolved signal to line up with spike-times
                    if numel(shifty) < numel(bins)
                        % pad with zeros
                        diffbins = numel(bins) - numel(shifty);
                        shifty = [zeros(diffbins,1); shifty];
                    end % can occasionally happen with width pars that are not integer multiples of step-size
                    spkfcn(:,j) = shifty;
                catch
                    keyboard
                end
                % keyboard
            case 'exponential'
                y = conv(h,spk);
                spkfcn(:,j) = y(1:end-shiftbase+1);   % shift convolved signal to line up with spike-times
        end
    end
    % keyboard
end

try
    % if not firing, then strip out cells
    spkfcn(:,nspikes==0) = [];
    idxs = Didxs; idxs(nspikes==0) = [];
catch
    keyboard
end


function [Sxy] = constructS(spkfcn,Nidxs,opts)

switch opts.Dmeth
    case {'corr','corrcoef'}
        Sxy = eval([opts.Dmeth '(spkfcn);']);     % compute correlation
        if any(isnan(Sxy))
            %keyboard
            y = find(isnan(Sxy));
            Sxy(y) = 0; % remove NaNs
        end
        Sxy(Sxy < 0) = 0;   % rectify correlation coefficient = similarity....
        % place zeros on the diagonal: no self-connections allowed...
        Sxy(eye(Nidxs)==1) = 0;
    case 'cosine'
        try
            % Sxy = 1 - squareform(pdist(spkfcn','cosine')); % get angle between spike trains; places 1s on the diagonal
            Sxy = squareform(1-pdist(spkfcn','cosine')); % places 0s on the diagonal
        catch
            %% if distances very small, then pdist throws error
            warning('Cosine distance too small, choose different binsize')
            Sxy = zeros(Nidxs*(Nidxs-1)/2,1);    % Sxy is upper triangle..
            % keyboard
        end
        %Sxy = 1 - pdist(spkfcn{loop}','cosine'); % get angle between spike trains
        % keyboard
        % need to re-scale this as "dynamic range" very small - see
        % e.g. Fellous et al, 2004
        % Fellous et al 2004 dynamic range function
        %             for ti = 1:numel(opts.tau)
        %                     % spin round tau options, make histogram, choose tau
        %                     % with flattest histogram
        %                     SxyC = 1 ./ (1 + exp(-(Sxy - mean(Sxy))./opts.tau(ti))); % pick tau to "flatten" histogram of similarity...
        %                     hst = histc(SxyC,[0:1/50:1]);
        %                     stdhst(ti) = std(hst);
        %                     % keyboard
        %             end
        %             % choose tau with min std....
        %             minix = find(stdhst == min(stdhst));
        %             Sxy = 1 ./ (1 + exp(-(Sxy - mean(Sxy))./opts.tau(minix(1))));
        
   % case 'xcorr'
    otherwise
        error('Unknown pairwise similarity metric specified')
end
