function [chr,fileprefix,channel] = tz_parseyeastname(yeastname)
%TZ_PARSEYEASTNAME parse the combined yeast name.
%   CHR = TZ_PARSEYEASTNAME(YEASTNAME) parse the string returned from
%   TZ_YEATIMG2FEATNAME. It returns the chormosome name.
%   
%   [CHR,FILEPREFIX,CHANNEL] = TZ_PARSEYEASTNAME(...) also returns file 
%   prefix and channel of the [yeast image]. 
%   
%   See also

%   16-Feb-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('Exactly 1 argument is required')
end

dotPosition = find(yeastname=='.');
if ~isempty(dotPosition)
    nameWithNoExtension = yeastname(1:dotPosition(1)-1);
else
    nameWithNoExtension = yeastname;
end
separator = '_';
tokens = tz_strtok(nameWithNoExtension,separator);

chr = tokens{1};
fileprefix = tz_cell2str(tokens(2:end-1),separator);
channel = tokens{end};