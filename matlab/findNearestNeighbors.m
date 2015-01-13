function nearest_nbrs = findNearestNeighbors(coords, warpCartCoords)
% coords should be a structure containing cartCoords, neighbors, and triangles
% THIS FUNCTION DOES NOT WORK!!!

cartCoords = coords.cartCoords;
N = size(cartCoords,2);
nbrs = [0:N-1; coords.neighbors];

wN = size(warpCartCoords,2);

totalNbrs = zeros(N,1);
for n = 1:N
	totalNbrs(n) = length(find(nbrs(:,n)~=-99));
end

nearest_nbrs = zeros(N,1);
for j = 1:wN
	
	closestNbr = j-1;
	while 1
		currNbrs = nbrs(1:totalNbrs(closestNbr+1),closestNbr+1);
		currNbrs+1
		nbrsCartCoords = cartCoords(:, currNbrs+1);
		
		projections = warpCartCoords(:,j)'*nbrsCartCoords;
		[maxVal,I]=max(projections);
		
		if closestNbr == currNbrs(I)
			break;
		else
			closestNbr = currNbrs(I);
		end
	end
	nearest_nbrs(j) = closestNbr;
end