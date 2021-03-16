function [objfeats,classlabel] = ...
    tz_combobjfeat2cell(combobj,combclass,combcellidx)
%TZ_COMBOBJFEAT2CELL Convert combined object feaures to MC objects.
%   OBJFEATS = TZ_COMBOBJFEAT2CELL(COMBOBJ,COMBCLASS,COMBCELLIDX) returns
%   a one-level cell array of object feature matrices, which is converted
%   from the matrix of all object features. COMBCLASS and COMBCELLIDX are
%   the vectors of class labels and cell indices.
%   
%   [OBJFEATS,CLASSLABEL] = TZ_COMBOBJFEAT2CELL(...) also returns the one
%   level cell array of cell index vectors.
%
%   See also TZ_MCF2COMBOBJFEATS

%   09-MAR-2004 Initial write T. Zhao
%   31-OCT-2004 Modified T. Zhao
%       - add comments
%       - change function name tz_combfeat2cell-->tz_combobjfeat2cell
%   Copyright (c) Murphy Lab, Carnegie Mellon University


if nargin < 3
    error('Exactly 3 arguments are required')
end

caclass=ml_findclass(combclass);

nclass=length(caclass);
k=1;

for i=1:nclass
    clabel=caclass{i}(1);
    cacell=ml_findclass(combcellidx(combclass==clabel));
    ncell=length(cacell);
    
    for j=1:ncell
        objfeats{k}=[combobj(combclass==clabel & combcellidx==cacell{j}(1),:)];
        classlabel(k)=clabel;
        k=k+1;
    end
end

classlabel=classlabel';