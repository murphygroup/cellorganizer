function s = tz_addfield(s,field,v)
%TZ_ADDFIELD Add structure field.
%   S = TZ_ADDFIELD(S,FIELD,V) add a filed with the name FIELD and value
%   V. S could be an empty structure.
%   
%   See also RMFIELD

%   18-Nov-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 3
    error('Exactly 3 arguments are required')
end

if isempty(s)
    clear s;
end

eval(['s.' field '=v;']);