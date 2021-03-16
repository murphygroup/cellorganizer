function features = tz_combinefeats_df(features1,features2)
%TZ_COMBINEFEATS_DF Combine drug image features.
%   FEATURESS = TZ_COMBINEFEATS_DF(FEATURES1,FEATURES2) combine two drug
%   image features along columns.

%   14-JUL-2003  Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

nprot1 = size(features1,3);
ndrug1 = size(features1,2);

nprot2 = size(features2,3);
ndrug2 = size(features2,2);

if (nprot1~=nprot2) | (ndrug1~=ndrug2)
    warning('Size unmatched!');
    features = nan;
    return
end

nprot=nprot1;
ndrug=ndrug1;

for i=1:nprot
    for j=1:ndrug
        f1 = features1{1,j,i};
        f2 = features2{1,j,i};
        features{1,j,i}=tz_combinefeats(f1,f2);
        
        f1 = features1{2,j,i};
        f2 = features2{2,j,i};
        features{2,j,i}=tz_combinefeats(f1,f2);
    end
end
