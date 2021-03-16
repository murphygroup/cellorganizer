function y = tz_logitrnd(beta,x,n)
%TZ_LOGITRND Generate samples from a logistic model.
%   Y = TZ_LOGITRND(BETA,X) returns a vector of samples that is for the
%   model:
%       P(Y=1|X) = 1/(1+exp(-X*BETA))
%   
%   X is a MxN [feature matrix] and beta is an (N+1)x1 vector.
%
%   Y = TZ_LOGITRND(BETA,X,N) generates samples based on binomial
%   distribtions P(Y|X)~binomial(N,1/(1+exp(-X*BETA))).
%
%   See also

%   20-Jan-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 2
    error('2 or 3 arguments are required')
end

if ~exist('n','var')
    n = 1;
end

% ps = 1./(1+exp(-x*beta));
ps = tz_evallogistic(x,beta);
y = binornd(n,ps);