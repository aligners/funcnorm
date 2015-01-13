function [numNodes, coords] = parseSurfaceFile(surfFile, opsDataType, saveDataType, maxNbrs)
% FUNCTION [numNodes, coords] = parseSurfaceFile(surfFile [, opsDataType='single' [, saveDataType='single']])
%
% *** INPUT ARGUMENTS ***
%	surfFile:  The filename of the surface file
%		Usually located in the output directory (such as [outDir, 'standard', num2str(mmRes), 'mm_', fs_surf, '.asc'];
%			Where: 	mmRes is the resolution of the surface mesh in mm
%					fs_surf is freesurfer surface type (e.g., 'sphere.reg')
%	opsDataType:  The datatype to load coords as (e.g., 'single')
%	saveDataType:  The datatype to save data out to (e.g., 'single')
%
% *** OUTPUT ARGUMENTS ***
%	numNodes:  the number of nodes in the mesh
%	coords:  a coordinates structure containing the following elements:
%		.cartCoords:  a 3 x numNodes matrix of cartesian coordinates
%		.triangles:  a 3 x numTriangles matrix of triplets of node numbers
%		.neighbors:  a 6 x numNodes matrix of node numbers (specifying the 6 nearest-neighbors of each node in the mesh)
%
%
% This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
% Please see AUTHORS and LICENSE file in the project root directory

if nargin == 1
    opsDataType = 'single';
    saveDataType = 'single';
end

if ~exist('maxNbrs'); maxNbrs = 6; end;

fp = fopen(surfFile, 'r');
if fp == -1
	error(['Could not open surface file:  ', surfFile]);
end

% Ignore first line of the surface file
fgetl(fp);
l = fscanf(fp, '%d %d\n', 2);
numNodes = l(1);
numTris = l(2);

coords.cartCoords = cast(fscanf(fp, '%f %f %f 0\n', [3, numNodes]), opsDataType);
coords.triangles = cast(fscanf(fp, '%d %d %d 0\n', [3, numTris]), 'int32');

fclose(fp);

% Generate the neighbors
[coords.neighbors, nbrCounts] = generateNeighborsMatrixFromTriangles(coords.triangles, numNodes, maxNbrs);