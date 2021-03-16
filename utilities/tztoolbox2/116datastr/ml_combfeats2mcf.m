function features=ml_combfeats2mcf(combfeats,combclass)
%ML_COMBFEATS2MCF Change feature structure from combined to MCF.
%   FEATURES = ML_COMBFEATS2MCF(COMBFEATS,COMBCLASS) returns MCF from the
%   combined feature matrix COMBFEATS and class labels COMBCLASS, which
%   is a column vector. COMBFEATS and COMBCLASS must have the same number
%   of rows.
%   
%   See also ML_MCF2COMBFEATS

%   ??-???-???? Initial write T. Zhao
%   31-OCT-2004 Modified T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

features=ml_findclass([combfeats,combclass]);
nclass=length(features);

for i=1:nclass
    features{i}=features{i}(:,1:end-2);
end