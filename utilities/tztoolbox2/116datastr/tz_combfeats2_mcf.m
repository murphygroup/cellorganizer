function feats=tz_combfeats2_mcf(feats1,feats2)
%TZ_COMBFEATS2_MCF Combine 2 MCFs. (Obsolete)
%
%See also TZ_COMBMCF2

%   FEATS = TZ_CLASSOBJCOM(FEATS1,FEATS2) returns a cell array of feature
%   matrices. Each element of the cell array is the combination of the
%   corresponding elements in FEATS1 and FEATS2 along columns. This is
%   supposed to combine different feature set for the same data set.

%   ??-???-???? Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

error(tz_genmsg('of','tz_combfeats2_mcf','tz_combmcf2');

if nargin < 2
    error('Exactly 2 arguments are required')
end

nclass = length(feats1);

for i=1:nclass
    feats{i}=[feats1{i},feats2{i}];
end
