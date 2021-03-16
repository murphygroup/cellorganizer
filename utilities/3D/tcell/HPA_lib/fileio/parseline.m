function output = parseline( strline, delimiter, negator)

if ~exist('delimiter','var')
    delimiter = ',';
end
if ~exist('negator','var')
    negator = '"';
end

idx = findstr(strline, delimiter);

idx2 = findstr(strline, negator);
ind = [];
for i=1:2:length(idx2)
    ind = [ind find((idx>idx2(i)) & (idx<idx2(i+1)))];
end
idx(ind) = [];


start = 1;
output = [];
for i=1:length(idx)
    output{i} = strline(start:idx(i)-1);
    start = idx(i)+1;
end

if max(idx)<length(strline)
    output{end+1} = strline( idx(end)+1:end);
end

return
