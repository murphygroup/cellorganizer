function allfeats=tz_cellnormdistfeat_mcf(alldists,allweights,s)
%TZ_CELLNORMDISTFEAT_MCF Calculate distance features for mcf cells.
%   ALLFEATS = TZ_CELLNORMDISTFEAT_MCF(ALLDISTS,ALLWEIGHTS,S) returns a
%   cell array of feature matrices based on the normalized distances
%   ALLDISTS and object weights ALLWEIGHTS. S specifies the structure of
%   ALLDISTS. If S is 1, ALLDISTS has a structure like 
%   {class}{cell}{[dist],...[dist]}, otherwise it is like
%   {class}{cell}[dist]. It is a bit different with the function 
%   TZ_CALCMORPHDISTFEATS_MCF, which does not consider weights and s .
%
%   ALLFEATS = TZ_CELLNORMDISTFEAT_MCF(ALLDISTS,ALLWEIGHTS) uses the
%   default value 1 of S.

%   ??-???-???? Initial write T. Zhao
%   31-OCT-2004 Modified T. Zhao
%       -a dd comments
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('2 or 3 arguments are required')
end

if ~exist('s','var')
    s=1;
end

nclass=length(alldists);

for i=1:nclass
    ncell=length(alldists{i});
    feats=[];
    for j=1:ncell
        nobj=length(alldists{i}{j});
        normdists=[];
        weights=[];
        for k=1:nobj
            if s==1
                normdists=[normdists;alldists{i}{j}{k}];
            else
                normdists=[normdists;alldists{i}{j}(k)];
            end
            if ~isempty(allweights)
                weights=[weights;allweights{i}{j}{k}(:,end)];
            end
        end
        feats(j,:)=tz_calcmorphdistfeats(normdists,weights);
    end
    allfeats{i}=feats;
end