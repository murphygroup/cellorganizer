function [combobjfeats,combclass,combcellidx] = ...
    tz_objfeatcomb2mcf(objfeats)
%TZ_OBJFEATCOMB2MCF Convert MCF object features to combined features.
%   COMBOBJFEATS = TZ_OBJFEATCOMB2MCF(OBJFEATS) returns the whole feature
%   matrix of all of the object features, which is the combination of the
%   two-level cell array of object features.
%   
%   [COMBOBJFEATS,COMBCLASS,COMBCELLIDX] = TZ_OBJFEATCOMB2MCF(...) also
%   returns combined class labels and combined cell indices.

%   15-Sep-2005 Initial write T. Zhao
%   ??-???-???? Initial write T. Zhao
%   03-NOV-2004 Modified T. Zhao
%       - add comments
%       - change function name tz_combineobjfeats_mcf --> tz_objfeatcomb2mcf
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argument is required')
end

nclass=length(objfeats);
combobjfeats=[];
combclass=[];
combcellidx=[];

for i=1:nclass
    ncell=length(objfeats{i});
    for j=1:ncell
        nobj=size(objfeats{i}{j},1);
        combobjfeats=[combobjfeats;objfeats{i}{j}];
        combclass=[combclass;zeros(nobj,1)+i];
        combcellidx=[combcellidx;zeros(nobj,1)+j];
    end
end