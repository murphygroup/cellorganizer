function indices = tz_findmulti(x,y)
%TZ_FINDMULTI Find a vector from another vector.
%   INDICES = TZ_FINDMULTI(X,Y) returns a cell array of vectors of indices of
%   the occurrence of X in Y.
%   
%   See also

%   29-May-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 2
    error('Exactly 2 arguments are required')
end

for i=1:length(x)
    indices{i} = find(y==x(i));
end
