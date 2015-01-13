function [warpDistances, distP1, distP2] = script_computeWarpConsistency(lambda, passNum, zeroCorrected)

subjects = {'cb', 'dm', 'hj', 'kd', 'kl', 'mh', 'ph', 'rb', 'se', 'sm'};
nsubj = length(subjects);
rho = 100;

if zeroCorrected == 1
	strCorr = '_corrected';
else
	strCorr = '';
end

warpDirP1 = ['../surf_data/alignment_raidersP1_', num2str(lambda), '/tmp/'];
warpDirP2 = ['../surf_data/alignment_raidersP2_', num2str(lambda), '/tmp/'];

surfFile = sprintf('../surf_data/standard2mm_sphere.reg.asc');
[numNodesPerHem, coords] = parseSurfaceFile(surfFile);
numNodes = 2*numNodesPerHem;
coords.cartCoords = normalizeCartesianCoordinates(coords.cartCoords);
coords = generateMultiHemisphereCoordsDataset(coords, 2);

minMM = 6;
warpDistances = [];
distP1 = [];
distP2 = [];
warpP1 = zeros(3, numNodes);
warpP2 = zeros(3, numNodes);
for subjNum = 1:nsubj
	subj = subjects{subjNum};
	
	fp = fopen([warpDirP1, 'warp_', subj, '_lh_', num2str(passNum), strCorr, '.bin'], 'rb');
	warpP1(1:3, 1:numNodesPerHem) = fread(fp, [3, numNodesPerHem], 'float32');fclose(fp);
	fp = fopen([warpDirP1, 'warp_', subj, '_rh_', num2str(passNum), strCorr, '.bin'], 'rb');
	warpP1(1:3, (numNodesPerHem+1):numNodes) = fread(fp, [3, numNodesPerHem], 'float32');fclose(fp);
	
	fp = fopen([warpDirP2, 'warp_', subj, '_lh_', num2str(passNum), strCorr, '.bin'], 'rb');
	warpP2(1:3, 1:numNodesPerHem) = fread(fp, [3, numNodesPerHem], 'float32');fclose(fp);
	fp = fopen([warpDirP2, 'warp_', subj, '_rh_', num2str(passNum), strCorr, '.bin'], 'rb');
	warpP2(1:3, (numNodesPerHem+1):numNodes) = fread(fp, [3, numNodesPerHem], 'float32');fclose(fp);

	coordsP1 = coords.cartCoords + warpP1;
	coordsP2 = coords.cartCoords + warpP2;
	
	d1 = computeGeodesicDistances(coords.cartCoords, coordsP1, rho);
	d2 = computeGeodesicDistances(coords.cartCoords, coordsP2, rho);
	locs = intersect(find(d1>minMM), find(d2>minMM));
	distP1 = [distP1, d1(locs)];
	distP2 = [distP2, d2(locs)];
	warpDistances = [warpDistances, computeGeodesicDistances(coordsP1(:, locs), coordsP2(:, locs), rho)];
end
