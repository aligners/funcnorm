function [warpDistsHist, xout] = script_computeWithinSubjWarpDistances(lambda, mode)

subjects = {'cb', 'dm', 'hj', 'kd', 'kl', 'mh', 'ph', 'rb', 'se', 'sm'};
nsubj = length(subjects);

if strcmp(mode, 'acrossSubj')
	warpDir = ['../surf_data/alignment_acrossSubjReg_raidersP1_', num2str(lambda), '/tmp/'];
	maxMM = 8;
	nbins=200;
elseif strcmp(mode, 'withinSubj')
	warpDir = ['../surf_data/alignment_withinSubjReg_raidersP1_', num2str(lambda), '/tmp/'];
	maxMM = 8;
	nbins=200;
else
	error(['improper mode ', mode]);
end

surfFile = '../surf_data/standard2mm_sphere.reg.asc';
[numNodesPerHem, coords] = parseSurfaceFile(surfFile);
numNodes = 2*numNodesPerHem;
coords.cartCoords = normalizeCartesianCoordinates(coords.cartCoords);
coords = generateMultiHemisphereCoordsDataset(coords, 2);

rho = 100;
warpDistsHist = zeros(1, nbins);
numNodesIncluded = 0;
freqLessThan1 = 0;
freqLessThan2 = 0;
totalNodes = 0;
warpDistMean = 0;
warpDistSSQ = 0;
if strcmp(mode, 'acrossSubj')
	runNum = 1;
	while 1
		fileNameLH = sprintf('%s/warp_runNum_%d_lh.bin', warpDir, runNum);
		fileNameRH = sprintf('%s/warp_runNum_%d_rh.bin', warpDir, runNum);
		
		if (exist(fileNameLH) == 2) && (exist(fileNameRH) == 2)
			disp(sprintf('Found %d', runNum));pause(0.001);
			fp = fopen(fileNameLH, 'rb');
			warp(1:3, 1:numNodesPerHem) = fread(fp, [3, numNodesPerHem], 'float32');fclose(fp);
			fp = fopen(fileNameRH, 'rb');
			warp(1:3, (numNodesPerHem+1):numNodes) = fread(fp, [3, numNodesPerHem], 'float32');fclose(fp);
		
			warpDists = computeGeodesicDistances(coords.cartCoords, coords.cartCoords+warp, rho);

            freqLessThan1 = freqLessThan1 + length(find(warpDists<1));
            freqLessThan2 = freqLessThan2 + length(find(warpDists<2));
            totalNodes = totalNodes + length(warpDists);

            locs = find(warpDists < maxMM);
			warpDistMean = warpDistMean + sum(warpDists(locs));
			warpDistSSQ = warpDistSSQ + sum(warpDists(locs).^2);
			numNodesIncluded = numNodesIncluded + length(locs);
			[w, xout] = hist(warpDists(locs), nbins);
			warpDistsHist = warpDistsHist + w;
			
			runNum = runNum + 1;
		else
			break;
		end
	end
else
	for subjNum = 1:nsubj
		subj = subjects{subjNum};
	
		suffix = '1stto2nd';
		fileNameLH = sprintf('%s/warp_%s_lh_%s.bin', warpDir, subj, suffix);
		fileNameRH = sprintf('%s/warp_%s_rh_%s.bin', warpDir, subj, suffix);
		if (exist(fileNameLH) == 2) && (exist(fileNameRH) == 2)
			disp(sprintf('Found %s %s', subj, suffix));pause(0.001);
			fp = fopen(fileNameLH, 'rb');
			warp(1:3, 1:numNodesPerHem) = fread(fp, [3, numNodesPerHem], 'float32');fclose(fp);
			fp = fopen(fileNameRH, 'rb');
			warp(1:3, (numNodesPerHem+1):numNodes) = fread(fp, [3, numNodesPerHem], 'float32');fclose(fp);
		
			warpDists = computeGeodesicDistances(coords.cartCoords, coords.cartCoords+warp, rho);

            freqLessThan1 = freqLessThan1 + length(find(warpDists<1));
            freqLessThan2 = freqLessThan2 + length(find(warpDists<2));
            totalNodes = totalNodes + length(warpDists);

            locs = find(warpDists < maxMM);
			warpDistMean = warpDistMean + sum(warpDists(locs));
			warpDistSSQ = warpDistSSQ + sum(warpDists(locs).^2);
			numNodesIncluded = numNodesIncluded + length(locs);
			[w, xout] = hist(warpDists(locs), nbins);
			warpDistsHist = warpDistsHist + w;
		end


		suffix = '2ndto1st';
		fileNameLH = sprintf('%s/warp_%s_lh_%s.bin', warpDir, subj, suffix);
		fileNameRH = sprintf('%s/warp_%s_rh_%s.bin', warpDir, subj, suffix);
		if (exist(fileNameLH) == 2) && (exist(fileNameRH) == 2)
			disp(sprintf('Found %s %s', subj, suffix));pause(0.001);
			fp = fopen(fileNameLH, 'rb');
			warp(1:3, 1:numNodesPerHem) = fread(fp, [3, numNodesPerHem], 'float32');fclose(fp);
			fp = fopen(fileNameRH, 'rb');
			warp(1:3, (numNodesPerHem+1):numNodes) = fread(fp, [3, numNodesPerHem], 'float32');fclose(fp);
		
			warpDists = computeGeodesicDistances(coords.cartCoords, coords.cartCoords+warp, rho);

            freqLessThan1 = freqLessThan1 + length(find(warpDists<1));
            freqLessThan2 = freqLessThan2 + length(find(warpDists<2));
            totalNodes = totalNodes + length(warpDists);
            
            locs = find(warpDists < maxMM);
			warpDistMean = warpDistMean + sum(warpDists(locs));
			warpDistSSQ = warpDistSSQ + sum(warpDists(locs).^2);
			numNodesIncluded = numNodesIncluded + length(locs);
			[w, xout] = hist(warpDists(locs), nbins);
			warpDistsHist = warpDistsHist + w;
		end
	end
end

if numNodesIncluded > 0
	warpDistMean = warpDistMean/numNodesIncluded;
	warpDistStdev = sqrt(warpDistSSQ/numNodesIncluded - warpDistMean^2);
	disp(['Warp Distance Mean:  ', num2str(warpDistMean)]);
	disp(['Warp Dist stdev:  ', num2str(warpDistStdev)]);
	warpDistsHist = warpDistsHist/numNodesIncluded;
	bar(xout, warpDistsHist);
	xlabel('Warp Distance (mm)');
	ylabel('Frequency');
	
	
%	locs = find(xout <= 1);
%	freqLessThan1 = sum(warpDistsHist(locs));
    freqLessThan1 = freqLessThan1/totalNodes;
	disp(['Frequency of nodes less than 1mm:  ', num2str(freqLessThan1)]);pause(0.001);
	
%	locs = find(xout <= 2);
%	freqLessThan2 = sum(warpDistsHist(locs));
    freqLessThan2 = freqLessThan2/totalNodes;
	disp(['Frequency of nodes less than 2mm:  ', num2str(freqLessThan2)]);pause(0.001);
	
%	locs = find(xout > 2);
%	freqGreaterThan2 = sum(warpDistsHist(locs));
	disp(['Frequency of nodes greater than 2mm:  ', num2str(1-freqLessThan2)]);pause(0.001);
	
%	figure;
%	plot(xout, warpDistsHist);
else
	error('Didn''t find anything!');
end