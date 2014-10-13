function [grps,Vpos,B,Q] = allevsplit(A,varargin)

% ALLEVSPLIT partition graph using all positive eigenvectors of modularity matrix
%   [S,U,B,Q] = ALLEVSPLIT(A) splits the vertices of the graph in adjacency matrix
%   A into multiple groups, defined in the vector S (each element taking an integer 
%   value to indicate group membership of the vertex at that index). It
%   also returns U, the vectors of all positive eigenvectors of the modularity matrix; B
%   the modularity matrix; and Q the modularity score for each repeat of the k-means
%   clustering. Applies suggestion of Reichardt & Bornhaldt
%   (2006) to work with directed networks.
%
%   ...= ALLEVSPLIT(...,DIST,N) sets k-means distance metric to string DIST (set to ''
%   to omit - default is squared Euclidean distance); runs the k-means
%   clustering N times (default is 10);
%
%   Notes: 
%   (1) This is a one-step multiple partition method, following up a
%   suggestion in Newman (2006) that all eigenvectors corresponding to positive
%   eigenvalues of the modularity matrix contain information about group
%   structure. The algorithm implemented takes the C such eigenvectors, and
%   uses k-means to clusters the nodes into k = C+1 groups.
%   (2) The k-means clustering is currently repeated 10 times for
%   robustness; clustering with highest Q score is retained
%   (3) Uses k-means++ algorithm to find the initial centroids for each
%   k-means test
%
%   References: 
%   (1) Newman, M. E. J. (2006) "Finding community structure in
%   networks using the eigenvectors of matrices". Phys Rev E, 74, 036104.
%
%   (2) Reichardt & Bornhaldt (2006) "Statistical mechanics of community detection".
%   Phys Rev E. 74, 016110
%
%   Mark Humphries 23/11/2009

strDist = 'sqEuclidean';
nreps = 10;
if nargin >= 2
    if ~isempty(varargin{1}) strDist = varargin{1}; end  
    if ~isempty(varargin{2}) nreps = varargin{2}; end
end

[n c] = size(A);
Abar = (A + A') / 2;

% generate expected graph
P = expectedA(A);
Pbar = (P + P') / 2;

% create modularity matrix
B = Abar - Pbar;

% keyboard

% find eigenvalues (D) and eigenvectors (V) of B
try
    [V,D] = eig(B);         % assumes symmetry
catch
    %keyboard
     V = zeros(size(B));
     D = zeros(size(B));
end

eg = diag(D);
if ~isreal(eg)
    warning('Complex eigenvalues have occurred')
    ns = sum(V > 0);
end

% eigenvectors corresponding to positive eigenvalues
egpos = find(eg > 1e-3);    % allow for rounding error
ngrps = numel(egpos) + 1;   % upper bound is one more group than vectors
Vpos = V(:,egpos);

% keyboard

ndivs = numel(2:ngrps);
Qmax = zeros(ndivs,1); Q = zeros(nreps,ndivs); qD = zeros(nreps,ndivs);
allgrps = zeros(n,nreps); grps = zeros(n,ndivs); Dgrps = zeros(n,ndivs);
C = cell(nreps,1); sumD = cell(nreps,1); allD = cell(nreps,1); 
Dmax = zeros(ndivs,1);

% run over all possible numbers of groups to max, and record groups with
% max Q....
for j = 1:ndivs

    %keyboard
    for rep = 1:nreps
        % randn('state',rep); rand('state',rep)
        try
            % j+1 is the number of groups this time...
            cpos = kmeansplus(Vpos,j+1);
            
            % [allgrps(:,rep),C{rep},sumD{rep},allD{rep}] = kmeans(Vpos,j+1,'Distance',strDist);
            [allgrps(:,rep),C{rep},sumD{rep},allD{rep}] = kmeans(Vpos,j+1,'Distance',strDist,'Start',cpos);
            
            % construct S matrix of group membership: each column is a group
            % See: Newman (2006) Eq 31
            S = zeros(n,ngrps);

            for loop = 1:ngrps
                S(:,loop) = (allgrps(:,rep) == loop);
            end

            % compute modularity
            Q(rep,j) = trace(S' * B * S);
            
            %% Compute Fellous's D index as indicator of k-means quality...
%             for k = 1:j+1
%                 notthisgrp = 1:n;
%                 thisgrp = find(allgrps(:,rep) == k);
%                 notthisgrp(thisgrp) = [];   % complementary set
%                 Dk(k) = numel(thisgrp) / (n-numel(thisgrp)) * (sum(allD{rep}(notthisgrp,k)) ./ sum(allD{rep}(thisgrp,k))); 
%                 if isinf(Dk(k))
%                     % a single point in group, so has zero distance
%                     Dk(k) = 0; % probably a pathological grouping...
%                 end
%             end
%             qD(rep,j) = sum(Dk) / (j+1);

        catch
            % if kmeans throws a wobbly, set to "no groups"...
            allgrps(:,rep) = zeros(n,1);
            Q(rep,j) = -1;
            qD(rep,j) = -1;
            % keyboard
        end

        % return group structure from max Q
        if Q(rep,j) > Qmax(j) 
            Qmax(j) = Q(rep,j); 
            grps(:,j) = allgrps(:,rep); % this is the returned group....
        end
        
        %keyboard
        % if using D score....
        % if clustering of better "quality" then keep this grouping...         
%         if qD(rep,j) > Dmax(j)
%             Dmax(j) = qD(rep,j); 
%             Dgrps(:,j) = allgrps(:,rep); % this is the returned group....
%         end
        
    end % end k-means loop
    
    % if negative modularity, then no groups!
    if Qmax(j) <= 0
        grps(:,j) = zeros(n,1);
    end

    % safety-catch: if consecutive attempts at groupings find no groups 
    % after multiple k-means, then extremely unlikely that adding
    % increasing k further (recruiting smaller +ve eigenvalues) will
    % find anything
    if j > 2 
        if Qmax(j) <= 0 & Qmax(j-1) <= 0
            break
        end
    end
end % end groups loop

if ~any(Qmax) % i.e. all less than zero...
    grps = zeros(n,1);  % no groups
else
    % return based on Q
    maxix = find(Qmax == max(Qmax));
    grps = grps(:,maxix);
    
    % return based on D
%     maxix = find(Dmax == max(Dmax));
%     grps = Dgrps(:,maxix);
end

% keyboard

function cpos = kmeansplus(pos,k)
    % pos = position of each point; k = number of clusters
    [r c] = size(pos); if c>r pos = pos'; end
    [np c] = size(pos);
    
    cd = zeros(k,1); % array of center point indexes
    % np = numel(pos) / k % number of data-points
    I = randperm(np); 
    cd(1) = I(1);   % pick one uniformly at random 
    for loop = 1:k-1
        % find all distances to current center
        d = nthroot(sum((repmat(pos(cd(loop),:),np,1) - pos)'.^2),k)';
        pc = d.^2 / sum(d.^2);  % select next center with probability proportional to distance squared...
        pc(cd(cd~=0)) = 0;    % set all existing centers to probability zero...
        cd(loop+1) = discreteinvrnd(pc,1,1); % select next center  
    end
    cpos = pos(cd,:);


