function [n,range] = tz_countnum(x)
%TZ_COUNTNUM Obsolete. See ML_COUNTNUM.
%   N = TZ_COUNTNUM(X) returns a vector which contains the number of integers
%   occurring in the vector X, which contains integers. n(1) is the number of
%   the minmum in X and n(end) is the number of maxmum in X. 
%   
%   [N,RANGE] = TZ_COUNTNUM(X) also returns the range of X.
%    
%   See also

%   28-May-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

error(tz_genmsg('of','tz_countnum','ml_countnum'));

if nargin < 1
    error('Exactly 1 argument is required')
end

range(1) = min(x);
range(2) = max(x);

n = hist(x,range(2)-range(1)+1);
