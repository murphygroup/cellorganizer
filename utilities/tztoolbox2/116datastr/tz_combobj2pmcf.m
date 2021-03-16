function [objfeats,classlabel] = ...
    tz_combobj2pmcf(combobj,combclass,combcellidx,labeled)
%TZ_COMBOBJ2MPCF Reorganize combined object features.
%   OBJFEATS = TZ_COMBOBJ2MPCF(COMBOBJ,COMBCLASS,COMBCELLIDX) returns 
%   one-level cell array of the object features. So it is called
%   PMCF which stands for partial MCF.
%   
%   OBJFEATS = TZ_COMBOBJ2MPCF(COMBOBJ,COMBCLASS,COMBCELLIDX,LABELED)
%   returns the cell array of the labeled feature matrices, in which
%   the last columns are cell indices.
%   
%   [OBJFEATS,CLASSLABEL] = TZ_COMBOBJ2MPCF(...) also returns the vector
%   of the sorted class labels

%   09-MAR-2004 Initial write T. Zhao
%   05-NOV-2004 Modified T. Zhao
%       - change function name tz_reorgobj --> tz_combobj2pmcf
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 3
    error('3 or 4 arguments are required')
end

if ~exist('labeled','var')
    labeled=1;
end

caclass=tz_findclass(combclass);

nclass=length(caclass);

for i=1:nclass
    classlabel(i)=caclass{i}(1);
    if labeled==1
        objfeats{i}=[combobj(find(combclass==classlabel(i)),:),combcellidx(find(combclass==classlabel(i)))];
    else
        objfeats{i}=[combobj(find(combclass==classlabel(i)),:)];
    end
end

classlabel=classlabel';

