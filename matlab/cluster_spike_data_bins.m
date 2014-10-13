function [G,grpsizes,ngrps,Beta,varargout] = cluster_spike_data_bins(spkdata,Didxs,binsize,T,varargin)

% CLUSTER_SPIKE_DATA_BINS cluster spike-data using time-bin based metrics
% [G,GS,N,C] = CLUSTER_SPIKE_DATA_BINS(S,I,B,T) clusters spike trains into
% groups, based on similarities between fixed time-width binned versions of
% the spike-trains. Similarities between all pairs are stored in similarity matrix, which forms
% the basis for the clustering.
%
% S is a 2-column vector: first column is ID stamps; second column is
% time-stamps. IDs can be either for events (e.g. stimulus presentations
% for single unit data; time-stamps are then relative to stimulus presentation) 
% or neurons (for multi-unit data; time-stamps are then relative to recording period).
%
% I is a column vector of all possible ID stamps - note that they need not all appear in S
% because the ID event may have caused no spikes or been a silent neuron!
%
% B is vector of binsizes to try; 
%
% T is 2-element vector of [start end] for all time-stamps (e.g. if using zero as 
% time-stamp for stimulus presentation: for 2 seconds of data starting at stimulus presentation T = [0 2]; for 2
% seconds starting 0.5 second before stimulus T = [-0.5 1.5];
%
% ... = CLUSTER_SPIKE_DATA_BINS(...,OPTS) sets all options contained in the
% OPTS structure - omitting a field will use the default option, indicated
% by the parentheses:
%       OPTS.Dmeth = {('Hamming') | 'cosine' | 'rCxy' | 'corrcoef' | 'spearman'} 
%           sets the similarity method to: normalised Hamming error (as used in Humphries, 2010); the cosine 
%           of angle between vectors; or rectifed Pearson's correlation coefficient
%           (so that is falls in [0,1]). Also included are
%           options to invoke correlation measures: linear
%           (Pearson's) correlation coefficient; Spearman's rank coefficient (no good for sparse data - i.e.
%           small bins). Both of these are for future development of the
%           algorithm with a different null model.
%       OPTS.Bln = {(0)|1} if set to 1, will convert the binned
%           spike-counts into binary arrays [choosing Hamming above does this by anyway]
%       OPTS.blnZscore = {(0)|1} if set to 1, will convert the binned spike
%           counts into Z-scores (number of standard deviations from mean). Note: if selected
%           will automatically override the binary vector representation if
%           both are selected. Is ignored if Hamming is chosen as
%           similarity
%       OPTS.blnGraph = {(0)|1} : converts similarity matrix into adjacency
%           matrix using permutation test to find p-value for similarity
%           between each pair of spike-trains
%       OPTS.alpha (0.05) : sets the significance threshold for graph conversion (has no effect otherwise)
%           pairs of spike-trains with p < alpha are considered significantly similar, 
%           so are assigned a 1 in the adjacency matrix, everything else is assigned a 0.       
%       OPTS.npermutes (100) : number of permutations (shuffles of each
%           train) to get p-values; if sufficient computational power is
%           available, then 1000 is a good starting point
%       OPTS.Qt (0.001) : quantisation step (in seconds) of shuffled spike-trains (default is 1 ms)
%       OPTS.nlimit (6) : the minimum number of retained entries from graph conversion required to
%           proceed with the clustering analysis (default is 6: two groups
%           of 3 nodes each) [has no effect if blnGraph = 0].
%       OPTS.blnAll = {0|(1)} - uses the all positive eigenvector method to cluster the spike trains 
%           (function ALLEVSPLIT), as used in Humphries, 2010; this works
%           best on full comparison matrices (i.e. blnGraph=0); setting blnAll = 0 uses the leading
%           eigenvector method instead (function MULTILEADEVSPLIT), which
%           works better on adjacency matrices.
%       OPTS.modopts ({}) - passes anything in this cell array to the corresponding optional
%           arguments of the clustering functions (ALLEVSPLIT or MULTILEADEVSPLIT) - see
%           their help for details. For example, for ALLEVSPLIT, the first cell will be
%           passed into the first optional argument (a string), 
%           the second cell into the second optional argument (a number) 
%       OPTS.overlap (0): fraction of binsize by which bins overlap; default is for non-overlapping bins, 
%           maximum must be less than 1. A value of 0.5 will overlap bins by half their width.    
%       OPTS.blnS = {(0)|1} - if set to 1, does not do the clustering algorithm, but only computes 
%           the similarity matrices for each binsize: useful for just creating a set of such matrices for further processing.
%           [Default is 0]
%
% Returns: cell array of all group structure G, one cell per binsize
% tested; each with a two-column vector: first column is ID stamp, second column
% is group membership (integer value) - note that there is no guarantee all IDs will be assigned a group. 
% A cell array of group sizes GS; An array of the number of groups for each binsize, N; 
% and measure of grouping quality C for each binsize: for the ALLEVSPLIT algorithm (Humphries, 2010), 
% this just returns maximum modularity Q; for MULTIEVSPLIT, this returns metric in (Humphries et al 2009).
%
% [...,Sxy,BD,Axy] = CLUSTER_SPIKE_DATA_BINS(...) are optional outputs useful for further 
% processing: Sxy is a cell array of the similarity matrices; BD is a cell array of 
% matrices containing the binned spike-train vectors - each column is one 
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
% (2) The novel algorithm described in Humphries (2011) corresponds to the
% default choices in this function. Other options are included for further
% exploration of applications of this analysis technique.
%
% (3) A good starting point for binsize is the lowest of either the mean or
% median ISI taken over all the spike data: if the ISIs are symmetrically
% distributed, these will be the same; if left-skewed (lots of long ISIs)
% then mean will be lower than median; if right-skewed (lots of short
% ISIs), median will be lower than mean. Hence taking the lower of the two
% will guarantee the majority of bins will contain at most one spike.
%
% (4) Alternatively, try using a bin-size optimisation routine first, then
% running this function
%
% (5) If an ID stamp has no associated spikes (e.g. there was no response
% to the stimulation or an error in recording) then the train with that ID
% is omitted from the analysis: the groupings in G will not contain that ID
% stamp.
%
%
% References:
% (1) Humphries, M. D. (2010) Spike-train communities: finding groups of
% similar spike-trains. Submitted.
% (2) Humphries, Wood & Gurney (2009) Dopamine-modulated dynamic cell assemblies generated by the GABAergic striatal
% microcircuit, Neural Networks, 22, 1174-1188.  
% 
%
% Mark Humphries 12/06/2011

[r c] = size(Didxs);
if r==1 Didxs = Didxs'; end   % make column vector

% defaults  
opts.Dmeth = 'Hamming';     % use Hamming distance
opts.Bln = 0;        % don't use binary vector
opts.blnZscore = 0;  % don't Zscore transform the spike-counts   
opts.blnGraph = 0;   % don't make graph
opts.nlimit = 6;     % minimum number of retained nodes for graph...
opts.alpha = 0.05;    % signficance threshold for retaining in graph
opts.npermutes = 100;   % default number of permutations
opts.Qt = 0.001;     % 1 ms quantisation step for shuffled spike-trains   

opts.blnAll = 1;     % use all eigenvector method
opts.modopts = '';   % set no options for the all eigenvector method...
opts.overlap = 0;    % fraction of binwidths that overlap; default is 0
opts.blnS = 0;       % run full function, including clustering algorithm

if nargin >= 5
    if isstruct(varargin{1}) 
        tempopts = varargin{1}; 
        fnames = fieldnames(tempopts);
        for i = 1:length(fnames)
            opts = setfield(opts,fnames{i},getfield(tempopts,fnames{i}));
        end
    end
end

% storage
scount = cell(numel(binsize),1);
Sxy = cell(numel(binsize),1);
Rn = cell(numel(binsize),1);        % retained neurons...
grpsizes = cell(numel(binsize),1);
ngrps = zeros(numel(binsize,1));

%% analyse data
for loop = 1:numel(binsize)
    
    %% bin data-spike trains
    [scount{loop},idxs] = bin_spiketrains(spkdata,Didxs,binsize(loop),T,opts);
    Nidxs = numel(idxs);
    
    %% now compute selected similarity between those vectors 
    [Sxy{loop}] = constructS(scount{loop},Nidxs,opts);

        %% if constructing adjacency matrix, then run permutation tests and
    %% find graph structure...
    if opts.blnGraph
        % run permutation test
        for ctrl = 1:opts.npermutes
            % get shuffled trains
            ctrlspkts = shuffle_intervals(spkdata,T,opts.Qt);
            % bin control spike trains
            [ctrlscount,ctrlidxs] = bin_spiketrains(ctrlspkts,Didxs,binsize(loop),T,opts);

            %% now compute selected distance between those functions
            ctrlSxy{ctrl} = constructS(ctrlscount,Nidxs,opts);
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

    
    % do modularity on that graph
    % note that edges must count only existing/non-existing to detect the
    % component!
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
    
    % keyboard
    
    % remap from retained index count to ID stamps...
    G{loop} = [idxs(Rn{loop}) Sc];  % [ID stamp; Group membership]
    
    % compute "complexity" measure 
    nzSxy = Sxy{loop}(Sxy{loop}~=0);
        
    if opts.blnAll
        % then return max Q as the measure of grouping
        maxQ =  max(max(Q)); if isempty(maxQ) maxQ = 0; end
        Beta(loop) = maxQ;
    
    else
         %otherwise if used MULTILEADEVSPLIT, must use some other metric:
         %this uses the NN (2009) paper metric
         Beta(loop) = sum(grpsizes{loop})/Nidxs * ngrps(loop) * abs(min(nzSxy) - median(nzSxy)); 
    end
    
end % end loop over binsizes

varargout{1} = Sxy;
varargout{2} = scount;
varargout{3} = Acs;

function [scount,idxs] = bin_spiketrains(spkdata,Didxs,binsize,T,opts)
    
    Nidxs = numel(Didxs);
    
    if opts.overlap
         % make overlapping window times
        winstep = binsize * (1 - opts.overlap);
        winstrt = T(1):winstep:T(2)-binsize;
        nwins = numel(winstrt);
        scount = zeros(numel(nwins),Nidxs);
    else
        % discrete bins
        bins = T(1):binsize:T(2);   % construct vector of bins
        scount = zeros(numel(bins),Nidxs);
    end
    
    for j = 1:Nidxs
        currix = find(spkdata(:,1) == Didxs(j));
        nspikes(j) = numel(currix);
        if nspikes(j) > 0
            if opts.overlap
                % then loop over the overlapping bins, and sum spikes in each one
                for k = 1:nwins
                    winspks = find(spkdata(currix,2) > winstrt(k) & spkdata(currix,2) <= winstrt(k)+binsize);
                    scount(k,j) = numel(winspks);
                end
                % keyboard
            else
                % discrete bins edge-to-edge
                % use MATLAB's histc function for expediency
                scount(:,j) = histc(spkdata(currix,2),bins);
            end
        end
        if opts.blnZscore & isempty(strfind(opts.Dmeth,'Hamming'))
            scount(:,j) = (scount(:,j) - mean(scount(:,j))) ./ std(scount(:,j));
        end
    end
    
    % if not firing, then strip out cells
    scount(:,nspikes==0) = [];
    idxs = Didxs; idxs(nspikes==0) = [];
    
    
function [Sxy] = constructS(scount,Nidxs,opts)    

    % keyboard
    switch opts.Dmeth
        case 'Hamming'
            % disp('Hamming')
            % keyboard
            scount(scount > 1) = 1; % make vector binary (active / not-active)
            Sxy = squareform(1 - pdist(scount','hamming')); % get normalised Hamming distance (error)....  
            
        case {'rCxy'}
            if opts.Bln & ~opts.blnZscore % then turn into binary vector 
                scount(scount > 1) = 1; % make vector binary (active / not-active)
            end
            
            Sxy = corrcoef(scount);
            Sxy(Sxy<0) = 0; % rectify the correlation coefficient...
            Sxy(eye(Nidxs)==1) = 0;   % no self-connections allowed....

        case {'corr','corrcoef'}
            if opts.Bln & ~opts.blnZscore  % then turn into binary vector 
                scount(scount > 1) = 1; % make vector binary (active / not-active)
            end
            Sxy = eval([opts.Dmeth '(scount);']);   % just use normal Pearson's correlation coefficient - will need different null model P in ALLEVSPLIT
            Sxy (eye(Nidxs)==1) = 0;   % no self-connections allowed....

        case 'spearman'
            if opts.Bln & ~opts.blnZscore % then turn into binary vector 
                scount(scount > 1) = 1; % make vector binary (active / not-active)
            end
            % just use normal Spearman's rank - will need
            % different null model P in ALLEVSPLIT
            Sxy = squareform(1 - pdist(scount','spearman')); % get Spearman's rank coefficient (pdist returns 1-rank)
            
        case 'cosine'
            %%% no option to make binary vector here: doesn't make sense
            %%% for cosine similarity measure - may as well use Hamming...
            try 
                Sxy = squareform(1 - pdist(scount','cosine')); % get angle between spike trains  
            catch
                %% if distances very small, then pdist throws error
                warning('Cosine distance too small, choose larger binsize')
                Sxy  = zeros(Nidxs*(Nidxs-1)/2,1);    % Sxy is upper triangle..
                % keyboard
            end
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
            
        otherwise
            error('Do not recognise similarity metric');
    end


