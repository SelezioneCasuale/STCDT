function [grps,maxEV] = multileadevsplit(A,varargin)

% MULTILEADEVSPLIT partition graph using leading eigenvector of modularity matrix
%   [S,U] = MULTILEADEVSPLIT(A) splits the vertices of the graph in adjacency matrix
%   A into n groups, defined in the vector S (each element taking value
%   {1,...,n} to indicate group membership of the vertex at that index). It
%   also returns U, the leading eigenvector of the modularity matrix.
%
%   MULTILEADEVSPLIT(...,FLAG) where FLAG is a string featuring any
%   combination of:
%       'r': adds the refinement step after each sub-graph,
%               so that each split maximises Q.
%       'a': keeps all splits found, even if Q < 0
%
%
%   Notes: 
%   (1) this repeats the leading eigenvector split on each sub-graph
%   until the change in modularity is not positive...
%
%   (2) or could use singular value decomposition to find eigenvalues,
%   as standard MatLab routine EIG() does not converge for some sample
%   directed graphs (e.g. mac95) even though they are made symmetric.
%   (Note though that using SVD does not split the dolphin graph into 3)
%
%   References: 
%   (1) Newman, M. E. J. (2006a) "Finding community structure in
%       networks using the eigenvectors of matrices". Phys Rev E, in press.
%   (2) Newman, M. E. J. (2006b) "Modularity and community structure in
%       networks". PNAS, in press.
%   (3) Reichardt & Bornholdt (2006) "Statistical mechanics of community detection" 
%       Phys Rev E. 74, 016110
%   
%   Mark Humphries 27/1/2009

blnRefine = 0;
blnAll = 0;
if nargin >= 2 
    if strfind(varargin{1},'r') blnRefine = 1; end
    if strfind(varargin{1},'a') blnAll = 1; end
end

% apply correction to A for directed graphs. If graph undirected then this
% has no effect
Abar = (A + A') / 2;    

[n c] = size(A);
grps = zeros(n,1);
ks = sum(Abar);        
m = sum(ks)/2;         

[S,maxEV,B] = splitsubgraph(A,ks,m);

Q = S' * B * S;

% check modularity has increased
if ~blnAll & Q <= 1e-5    % very close or less than zero (deals with rounding errors)
    % if not increased then return nothing
    S = zeros(n,1);
    return
end

% make subgraphs
n1 = sum(S>0);
n2 = sum(S<=0);
subG1 = A(S>0,S>0);
subG2 = A(S<=0,S<=0);
grps = real(S>0) + 1;     % group 1 and group 2

subqueue = {subG1, subG2};
grpvertices = {find(S>0), find(S<=0)};

while ~isempty(subqueue)
    newA = subqueue{1};
    [S,maxEV,B] = splitsubgraph(newA,ks(grpvertices{1}),m); 
  
    if blnRefine
        S = refinesplit(S,B);
    end
      
    % compute Q
    % keyboard
    try
        Q = S' * B * S; 
    catch
        keyboard
    end
    
    %keyboard
    
    % split this group if:
    % (a) all splits are kept, and there is a split OR
    % (b) not all splits are kept, but Q > 0 (Newman 2006a suggestion)  
    if (blnAll & any(S>0) & any(S<0)) | Q > 0
        subG1 = newA(S>0,S>0);
        subG2 = newA(S<=0,S<=0);
        grp1 = grpvertices{1}(S>0);     % the ids of the original vertices
        grp2 = grpvertices{1}(S<=0);    % the ids of the original vertices
        grps(grp1) = max(grps)+1;       % re-number first group, second group retains numbers
        %grps(grp2) = max(grps)+1;
        % delete split graph from stack
        subqueue(1) = [];
        grpvertices(1) = [];
        % add new ones
        grpvertices = [{grp1,grp2} grpvertices];        
        subqueue = [{subG1,subG2} subqueue];
    else
        % don't split any further, delete top of stack
        subqueue(1) = [];
        grpvertices(1) = [];
    end

    % keyboard
end



function [S,maxEV,B] = splitsubgraph(A,ks,m) 
    % Note: if graph is directed, the correction to A should be
    % applied
    Abar = (A + A') / 2;   
    
    S = zeros(length(A),1);
    
    % generate expected graph
    P = expectedA(A);
    
    %keyboard
    
    % correct P
    Pbar = (P + P') / 2;
    
    kg = sum(Abar);        % in-degrees of sub-graph
    dg = sum(kg);       % total edges of sub-graph
    
    % this correction should apply without problem to a directed graph
    % because by this stage the graph will be undirected. The key is pass
    % the correct value of ks to this function.
    correction = kg - ks.*dg/(2*m);     % definitely 2*m here: is divided by 2 at start :-)
    
    % create modularity matrix - correction applied only to self-self
    % connections
    B = Abar - Pbar - diag(correction);   
        
    % find eigenvalues (D) and eigenvectors (V) of B
    [V,D] = eig(B);         % assumes symmetry
    eg = diag(D);
    if ~isreal(eg)
        warning('Complex eigenvalues have occurred')
    end
    
    % keyboard
    
     % largest positive eigenvalue is last element
    maxeig = max(eg); %% note - can have the same eigenvalue if nodes have no links in subgraph!
    ix = find(eg == maxeig(1));
    maxEV = V(:,ix(1)); % if more than one with the same max eigenvalue, just pick first...
    
        
    %%%% use singular value decomposition %%%%%%%%%%%%
%     [u,s,v] = svd(B);
%     sv = diag(s);
%     maxEV = v(:,find(sv == max(sv)));
     
    
    S(maxEV > 0) = 1;       % group 1
    S(maxEV <= 0) = -1;     % group 2                
    
    