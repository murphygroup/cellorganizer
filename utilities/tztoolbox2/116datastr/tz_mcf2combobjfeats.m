function [combfeats,combclass,combcellidx,combobjidx] = ...
    tz_mcf2combobjfeats(features)
%TZ_MCF2COMBOBJFEATS Conver mcf object featrures to combined one.
%   COMBFEATS = TZ_MCF2COMBOBJFEATS(FEATURES) returns the combined feature
%   of all of the object features in the two-level cell array.
%   
%   [COMBFEATS,COMBCLASS,COMBCELLIDX],COMBOBJIDX = TZ_MCF2COMBOBJFEATS(...)
%   also returns class labels, cell indices and object indices
%   
%   See also TZ_COMBOBJFEAT2CELL TZ_COMBOBJFEAT2MCF

%   13-NOV-2004 Initial write TING Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argument is required')
end

nclass=length(features);
combfeats=[];
combclass=[];
combcellidx=[];
combobjidx = [];

for i=1:nclass
    for j=1:length(features{i})
        combfeats=[combfeats;features{i}{j}];
        nobj=size(features{i}{j},1);
        combcellidx=[combcellidx;zeros(nobj,1)+j];
        combclass=[combclass;zeros(nobj,1)+i];
        combobjidx = [combobjidx;[1:nobj]'];
    end
end