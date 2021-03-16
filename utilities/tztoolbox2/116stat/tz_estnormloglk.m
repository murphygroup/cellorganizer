function lp = tz_estnormloglk(x,param)
%TZ_ESTNORMLOGLK Estimated log lilkelihood of data.
%   LP = TZ_ESTNORMLOGLK(X) returns the log likelihood of the [feature
%   matrix] X based on the estimated normal distribution of X.
%   
%   LP = TZ_ESTN ORMLOGLK(X,PARAM) let users specify how to estimate the
%   covariance matrix. PARAM is a structure with a field 'covmeth' to
%   specify covariance matrix estimation method:
%       'covmeth' - 
%           'mlk' : use matlab COV function
%           'ind' : assume all variables are independent
%
%   See also TZ_ESTNORMLOGLK

%   05-Jun-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('1 or 2 arguments are required')
end

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param,struct('covmeth','mlk'));

switch param.covmeth
    case 'mlk'
        s = cov(x);
    case 'ind'
        s = diag(var(x,0,1));
    otherwise
        error('unrecognized covariance estimation method');
end

mu = mean(x,1);

lp = tz_normloglk(x,mu,s);