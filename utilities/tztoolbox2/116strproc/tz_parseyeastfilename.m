function [orf,ver,x,exposure,filter,desginator] = ...
    tz_parseyeastfilename(filename)
%TZ_PARSEYEASTFILENAME Parse yeast file name.
%   ORF = TZ_PARSEYEASTFILENAME(FILENAME) returns open reading frame (ORF)
%   from the yeast file name FILENAME, which has the sturcture:
%       ORF_ver##-X-EXPOSURE_FILTER[_colocDESGINATOR
%   See http://yeastgfp.ucsf.edu for more details.   
%
%   [ORF,VER,X,EXPOSURE,FILTER,DESGINATOR] = TZ_PARSEYEASTFILENAME(...)
%   also returns VER,X,EXPOSURE,FILTER,DESGINATOR.
%   
%   See also

%   21-Feb-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('Exactly 1 argument is required')
end

filename = strtok(filename,'.');

tokens = tz_strtok(filename,'_');

orf = tokens{1};
tokens2 = tz_strtok(tokens{2},'-');

ver = tokens2{1};
x = tokens2{2};
exposure = tokens2{3};
filter = tokens{3};
if length(tokens)>3
    desginator = tokens{4};
else
    desginator = [];
end


