function lp = tz_normloglk(x,mu,s)
%TZ_NORMLOGLK Log likelihood of normal distribution.
%   LP = TZ_NORMLOGLK(X,MEAN,SIGMA) returns the log likelihood of the
%   [feature matrix] for multivariate normal distribtuion with mean MU and
%   covariance matrix S.
%   
%   See also

%   05-Jun-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 3
    error('Exactly 3 arguments are required')
end

nsample = size(x,1);
nvar = size(x,2);

lp = -nsample*nvar/2*log(2*pi)-nsample*log(det(s))/2;

x = ml_addrow(x,-mu)*sqrtm(inv(s));

lp = lp-sum(x(:).^2)/2;
