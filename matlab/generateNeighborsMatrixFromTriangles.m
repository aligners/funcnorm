function [nbrs, nbrCounts] = generateNeighborsMatrixFromTriangles(T, numNodes, maxNbrs)
% FUNCTION [nbrs, nbrCounts] = generateNeighborsMatrixFromTriangles(T, numNodes, maxNbrs)
% nbrs is maxNbrs x numNodes matrix
% nbrCounts is 1 x numNodes vector (number of neighbors of each node)
%
% T is 3 x numTriangles matrix
% numNodes is the total number of nodes
% maxNbrs is the maximum number of neighbors of any node

% This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
% Please see AUTHORS and LICENSE file in the project root directory

nbrs = -99*ones(maxNbrs, numNodes, 'int32');
nbrCounts = zeros(1, numNodes, 'int32');
numT = size(T, 2);

for j = 1:numT
    t1 = T(1, j);
    t2 = T(2, j);
    t3 = T(3, j);

    % First handle t1
    ind = t1 + 1;
    n1 = t2;
    n2 = t3;
    c = nbrCounts(ind);
    
    if length(find(nbrs(:, ind) == n1)) == 0
        c = c + 1;
        if c > maxNbrs
        	nbrs = [nbrs; -99*ones(1,numNodes, 'int32')];
        end
        nbrs(c, ind) = n1;
    end
    if length(find(nbrs(:, ind) == n2)) == 0
        c = c + 1;
        if c > maxNbrs
        	nbrs = [nbrs; -99*ones(1,numNodes, 'int32')];
        end
        nbrs(c, ind) = n2;
    end
    nbrCounts(ind) = c;

    
    % Now handle t2
    ind = t2 + 1;
    n1 = t1;
    n2 = t3;
    c = nbrCounts(ind);
    
    if length(find(nbrs(:, ind) == n1)) == 0
        c = c + 1;
        if c > maxNbrs
        	nbrs = [nbrs; -99*ones(1,numNodes, 'int32')];
        end
        nbrs(c, ind) = n1;
    end
    if length(find(nbrs(:, ind) == n2)) == 0
        c = c + 1;
        if c > maxNbrs
        	nbrs = [nbrs; -99*ones(1,numNodes, 'int32')];
        end
        nbrs(c, ind) = n2;
    end
    nbrCounts(ind) = c;

    % Now handle t3
    ind = t3 + 1;
    n1 = t1;
    n2 = t2;
    c = nbrCounts(ind);
    
    if length(find(nbrs(:, ind) == n1)) == 0
        c = c + 1;
        if c > maxNbrs
        	nbrs = [nbrs; -99*ones(1,numNodes, 'int32')];
        end
        nbrs(c, ind) = n1;
    end
    if length(find(nbrs(:, ind) == n2)) == 0
        c = c + 1;
        if c > maxNbrs
        	nbrs = [nbrs; -99*ones(1,numNodes, 'int32')];
        end
        nbrs(c, ind) = n2;
    end
    nbrCounts(ind) = c;    
end