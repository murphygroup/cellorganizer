function features=tz_combinemcf(features1,features2)
%TZ_COMBINEMCF Combine two MCFs.
%   FEATURESS = TZ_COMBINEMCF(FEATURES1,FEATURES2) combines two MCFs along
%   columns and returns the new MCF. Tt seems that it does the same thing 
%   as TZ_COMBFEATS2_MCF.

%   27-MAY-2004 Initial Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

if isempty(features1)
    features=features2;
    return;
end

if isempty(features2)
    features=features1;
    return;
end

if(length(features1)~=length(features2))
    error('class size not matched')
end

for i=1:length(features1)
    features{i}=tz_combinefeats(features1{i},features2{i});
end