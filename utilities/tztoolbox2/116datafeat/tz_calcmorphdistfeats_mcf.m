function allfeats=tz_calcmorphdistfeats_mcf(alldists)
%TZ_CALCMORPHDISTFEATS_MCF Calculate distance features for mcf.
%   ALLFEATS = TZ_CALCMORPHDISTFEATS_MCF(ALLDISTS) returns a cell array
%   of the feature matrices of a mcf cell array ALLDISTS.

%   ??-???-???? Initial write T. Zhao
%   30-OCT-2004 Modified T. Zhao
%       - add comments
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argument is required')
end

nclass=length(alldists);

for i=1:nclass
    ncell=length(alldists{i});
    feats=[];
    for j=1:ncell
        feats(j,:)=tz_calcmorphdistfeats(alldists{i}{j},[]);
    end
    allfeats{i}=feats;
end
    