function y = tz_wquantile(x,p)
%TZ_WQUANTILE Obsolete. See ML_WQUANTILE.
%   Y = TZ_WQUANTILE(X,P) returns the P quantile of the weights in X, which 
%   is a row vector or column vector.
%   
%   See also

%   26-Jan-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

error(tz_genmsg('of','tz_wquantile','ml_wquantile'));

if nargin < 2
    error('Exactly 2 arguments are required')
end

if size(x,1)>1
    x = x';
end

[s,idx] = sort(x,2,'descend');

cs = cumsum(s);
cs = cs/cs(end)>p;
[tmp,sidx] = max(cs);

y = s(sidx);






