function x = tz_gscross(mu1,sigma1,mu2,sigma2)
%TZ_GSCROSS Obsolete. See ML_GSCROSS.
%   X = TZ_GSCROSS(MU1,SIGMA1,MU2,SIGMA2) returns the values at which two
%   Gaussian distributions have equal density. The two Gaussian
%   distributions are:
%       X1 ~ N(MU1,SIGMA1)
%       X2 ~ N(MU2,SIGMA2)
%   
%   Document: normbound.doc. Ask T. Zhao if you need the document.
%
%   See also

%   27-Jan-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

error(tz_genmsg('of','tz_gscross','ml_gscross'));

if nargin < 4
    error('Exactly 4 arguments are required')
end

c(1) = (1/sigma2^2-1/sigma1^2)/2;
c(2) = -(mu2/sigma2^2-mu1/sigma1^2);
c(3) = (mu2^2/sigma2^2-mu1^2/sigma1^2)/2+log(sigma2/sigma1);

x = roots(c);
