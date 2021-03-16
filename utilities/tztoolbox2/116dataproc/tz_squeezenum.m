function k = tz_squeezenum(n)
%TZ_SQUEEZENUM Normalize integers.
%   K = TZ_SQUEEZENUM(N) returns a vector or matrix from N. K has the same
%   size as N but it only contains integers from 1 to M, where M is the
%   number of unique values in N.
%   
%   See also

%   18-Nov-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('Exactly 1 argument is required')
end

[a,b] = unique(n);
k=zeros(size(n));
for j=1:length(a)
    k(n==a(j)) = j;
end


