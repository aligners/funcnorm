function [data,num_nodes,num_timepoints] = read_niml_binary(filename, opsDataType, saveDataType)
%[data,num_nodes,num_timepoints] = read_niml_binary(filename, opsDataType, saveDataType)
% Reads in a niml binary file created by SUMA 3dVol2Surf -out_niml option
% warning: is kind of a hack resulting from examining only 2 such files!
%  opsDataType: defaults to 'single', is the precision of the output data
%  saveDataType:defaults to 'single', is the precision of the input binary
%                  (the header of the niml file should say "float")

% This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
% Please see AUTHORS and LICENSE file in the project root directory

if nargin == 1
    opsDataType = 'single';
    saveDataType = 'single';
end

fid=fopen(filename,'rb');
if (-1 == fid)
    error(['cannot open file:  ', filename]);
end

line=fgets(fid);
while isempty(findstr(line,'SPARSE_DATA'))
    line=fgets(fid);
end

line=fgetl(fid);
while isempty(findstr(line,'ni_form'))
    line=fgets(fid);
end

endian_loc=findstr(line,'sbfirst') - 1;

if line(endian_loc) == 'l'
    endian='l';
else
    endian='b';
end

line=fgets(fid);
while isempty(findstr(line,'ni_type'))
    line=fgets(fid);
end

quoteloc=findstr(line,'"');
asteriskloc=findstr(line,'*');

if isempty(asteriskloc)
    num_timepoints = 1;
else
    num_timepoints = str2num(line(quoteloc(1)+1:asteriskloc-1));
end

line=fgets(fid);
while isempty(findstr(line,'ni_dimen'))
    line=fgets(fid);
end

quoteloc=findstr(line,'"');
num_nodes = str2num(line(quoteloc(1)+1:quoteloc(2)-1));

fprev = ftell(fid);
line=fgets(fid);
loc=findstr(line,'>');
while isempty(loc)
    fprev = ftell(fid);
    line=fgets(fid);
    loc=findstr(line,'>');
end
if length(loc) > 1
    loc = loc(1);
end
fseek(fid, fprev+loc,'bof');

try
	data=fread(fid,[num_timepoints,num_nodes], [saveDataType, '=>', opsDataType], 0, endian);
catch exception
	filename
	[num_timepoints, num_nodes]
	throw(exception);
end

fclose(fid);

