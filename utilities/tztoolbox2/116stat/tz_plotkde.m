function tz_plotkde(x,varargin)
%TZ_PLOTKDE Plot estimated density from 1D data
%   TZ_PLOTKDE(X) plots the estimated density of X, which is a row or
%   column vector.
%   
%   See also

%   02-Nov-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('At least 1 argument is required')
end

if all(size(x)>1)
    error('X must be a vector');
end

if size(x,1)>1
    x = x';
end

mu = mean(x);
sigma = std(x);

y = tz_ppoints(mu-sigma*3,mu+sigma*3,100);

realpdf = normpdf(y,0,4);

trkde=kde(x,'lcv');
plot(y,evaluate(trkde,y),varargin{:});



