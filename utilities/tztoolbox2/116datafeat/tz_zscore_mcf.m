function normfeats=tz_zscore_mcf(features)
%TZ_ZSCORE_MCF Normalize MCF to mean 0 and variance 1.
%   NORMFEATS = TZ_ZSCORE_MCF(FEATURES) returns the normaized MCF by 
%   z-scoring, which is done on the whole data.

%   07-Aug-2004 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argument is required')
end

nclass=length(features);
combfeats=[];

for i=1:nclass
    combfeats=[combfeats;features{i}];
    ncell(i)=size(features{i},1);
end

normcombfeats=zscore(combfeats);

for i=1:nclass
    normfeats{i}=normcombfeats(1:ncell(i),:);
    normcombfeats(1:ncell(i),:)=[];
end
    
