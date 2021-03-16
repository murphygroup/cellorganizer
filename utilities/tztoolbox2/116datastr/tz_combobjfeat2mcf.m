function objfeats = tz_combobjfeat2mcf(combobjfeat,combclass,combcellidx)
%TZ_COMBOBJFEAT2MCF Convert an object feature matrix to a cell array.
%   OBJFEATS = TZ_COMBOBJFEAT2MCF(COMBOBJFEAT,COMBCLASS,COMBCELLIDX)
%   returns a 2-level cell array of [feature matrix]s. COMBOBJFEAT is a
%   [feature matrix], COMBCLASS is a [label vector] for class labels of the
%   features and COMBCELLIDX is a vector for cell indices in each class.
%   
%   See also TZ_COMBOBJFEAT2CELL TZ_MCF2COMBOBJFEATS

%   13-May-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 3
    error('Exactly 3 arguments are required')
end

nclass = max(combclass);
for i=1:nclass
    tmpfeat = combobjfeat(combclass==i,:);
    tmpcellidx = combcellidx(combclass==i);
    objfeats{i} = ml_combfeats2mcf(tmpfeat,tmpcellidx);
end