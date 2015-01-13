function num_nodes = getNumNodes(filename)

fid=fopen(filename,'rb');
if (-1 == fid)
    error('cannot open file');
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

num_timepoints = str2num(line(quoteloc(1)+1:asteriskloc-1));

line=fgets(fid);
while isempty(findstr(line,'ni_dimen'))
    line=fgets(fid);
end

quoteloc=findstr(line,'"');
num_nodes = str2num(line(quoteloc(1)+1:quoteloc(2)-1));
