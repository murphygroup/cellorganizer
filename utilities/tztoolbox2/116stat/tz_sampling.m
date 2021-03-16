function data = tz_sampling(f,n)
%TZ_SAMPLING Sampling data from a probability distribution.
%   DATA = TZ_SAMPLING(F,N) returns a vector or matrix which are drawn
%   from the distribution F, which is a structure (see TZ_TRAINLK). N is
%   the number of samples. It must be an integer. Each row of DATA is one
%   sample. F could also contain a field 'transform', which is a [general
%   function]. The data sampled from F will be transformed by the inverse
%   function of F.transform.
%   
%   See also TZ_TRAINLK

%   10-Jan-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 2
    error('Exactly 2 arguments are required')
end

switch f.method
    case 'pca'
        coeff = mvnrnd(f.coeffmean,f.coeffcov,n);
        data = ml_addrow(coeff*f.pcavec,f.mean);
    case 'mvn'
        data = mvnrnd(f.mu,f.sigma,n);
    case 'exp'
        data = exprnd(f.beta,n);
    case 'norm'
        data = normrnd(f.mu,f.sigma,n);

    otherwise
        error('Unrecognized probability funciton');
end

if isfield(f,'transform')
    data = ml_evalfun(data,ml_getinvfun(f.transform));
end