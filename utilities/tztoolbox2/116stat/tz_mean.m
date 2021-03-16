function mu = tz_mean(x,weights)
%TZ_MEAN Mean of data.
%   MU = TZ_MEAN(X) is the same as MEAN(X,1).
%   
%   MU = TZ_MEAN(X,WEIGHTS) uses weights for data when calculating the
%   mean.
%   
%   See also

%   17-Nov-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('1 or 2 arguments are required')
end

if ~exist('weights','var')
    mu = mean(x,1);
else
    weights = weights/sum(weights);
    mu=sum(repmat(weights,1,size(x,2)).*x,1);
end

