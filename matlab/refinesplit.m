function S = refinesplit(S,B,varargin)

% REFINESPLIT improve modularity maximisation
%   N = REFINESPLIT(S,B) attempts to further maximise the modularity value
%   of the network, described in modularity matrix B, given the current 2 group
%   partition in vector S. It returns the vector N of new groups that 
%   maximise modularity (values in {-1,+1}).
%
%   REFINESPLIT(...,'i') returns the groups in {1,2} instead of {-1,1},
%   which is more useful for post-processing.
%
%   Notes: 
%   #1: This works by taking each vertex in turn, assigning it to the
%   other group, and recomputing modularity Q. The group assignment that
%   maximises Q is kept. This is repeated until no further increases in Q
%   are obtained (with each vertex only moving once to a new group - i.e. 
%   once moved, they cannot move back!).
%
%   #2: Requires that group membership be +1 and -1, so converts the
%   membership vector S if required.
%
%   References: 
%   (1) Newman, M. E. J. (2006) "Finding community structure in
%       networks using the eigenvectors of matrices". Phys Rev E, in press.
%   (2) Newman, M. E. J. (2006) "Modularity and community structure in
%       networks". PNAS, in press.
%   
%   Mark Humphries 16/10/2006

% keyboard

Q = S' * B * S;
newQ = Q;
prevQ = Q;
blnSwap = 1;
G = sort(S);
g1 = -1;
g2 = 1;
if G(1) == G(end)
    % then no split to refine!
    return
end

S(S==G(1)) = g1;
S(S==G(end)) = g2;

[n c] = size(B);
moved = 0;

while blnSwap
    tempQ = zeros(n,1);
    for i = 1:n
        if i ~= moved
            Sp = S;
            if Sp(i) == g1 Sp(i) = g2;
            else Sp(i) = g1;
            end
            % recompute modularity
            tempQ(i) = Sp' * B * Sp; 
        else
            % compute modularity of unmoved vertex
            Sp = S;
            tempQ(i) = Sp' * B * Sp; 
        end
    end
    bestnode = find(tempQ == max(tempQ) & tempQ > newQ);
    if ~isempty(bestnode)
        % there could be more than one with exact same change in modularity
        % so take the first
        moved = [moved bestnode(1)];
    end
    
    if max(tempQ) > prevQ
        newQ = max(tempQ);
        prevQ = newQ;
        % swap permanently
        if S(bestnode) == g1 S(bestnode) = g2;
        else S(bestnode) = g1;
        end
    else 
        % no increase in Q
        blnSwap = 0;
    end
end

if nargin >= 3 & varargin{1} == 'i'
    S(S==g2) = 2;
    S(S==g1) = 1;
end

