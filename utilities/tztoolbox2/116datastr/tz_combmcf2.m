function feats=tz_combmcf2(feats1,feats2)
%TZ_COMBFEATS2 Combine 2 MCFs.
%   FEATS = TZ_COMBFEATS2(FEATS1,FEATS2) returns a cell array of feature
%   matrices. Each element of the cell array is the combination of the
%   corresponding elements in FEATS1 and FEATS2 along columns. This is
%   supposed to combine different feature set for the same data set.

%   ??-???-???? Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

nclass = length(feats1);

for i=1:nclass
    feats{i}=[feats1{i},feats2{i}];
end
