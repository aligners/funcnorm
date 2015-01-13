function num_timepoints = getNumTimepoints(filename)
% FUNCTION num_timepoints = getNumTimepoints(filename)
%
% *** INPUT ARGUMENTS ***
% 	filename (assumed to be .niml.dset)
%
% *** OUTPUT ARGUMENTS
%	num_timepoints:  the number of timepoints in the fMRI time-series
%
%
% This file is part of the Functional Normalization Toolbox, (c) 2008 by the authors.
% Please see AUTHORS and LICENSE file in the project root directory

fid=fopen(filename,'rb');
if (-1 == fid)
    error(['cannot open file:  ', filename]);
end

line=fgetl(fid);
while isempty(findstr(line,'SPARSE_DATA'))
    line=fgetl(fid);
end

line=fgetl(fid);
while isempty(findstr(line,'ni_form'))
    line=fgetl(fid);
end

endian_loc=findstr(line,'sbfirst') - 1;

if line(endian_loc) == 'l'
    endian='l';
else
    endian='b';
end

line=fgetl(fid);
while isempty(findstr(line,'ni_type'))
    line=fgetl(fid);
end

quoteloc=findstr(line,'"');
asteriskloc=findstr(line,'*');

num_timepoints = str2num(line(quoteloc(1)+1:asteriskloc-1));
