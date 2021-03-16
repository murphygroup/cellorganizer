function files2 = tz_sortfile(files,pos)
%TZ_SORTFILE Sort the file names with their numbers.
%   FILES2 = TZ_SORTFILE(FILES) sorts the string array FILES in ascending 
%   order. The orders are decided by the number in each string in FILES.
%
%   FILES2 = TZ_SORTFILE(FILES,POS) specifies the position wheret the 
%   number is obtained from each string.
%
%   See ML_GETFILENUM for more details.
%   
%   See also ML_GETFILENUM

%   20-Jan-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('1 or 2 arguments are required')
end

if ~exist('pos','var')
    pos = 1;
end

if isempty(files)
    files2 = {};
    return;
end

for i=1:length(files)
    fileNumbers(i) = ml_getfilenum(files{i},pos);
end

[tmp,sortIndices] = sort(fileNumbers);

files2 = files(sortIndices);
