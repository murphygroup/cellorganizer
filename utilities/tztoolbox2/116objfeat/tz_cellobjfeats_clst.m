function features = ...
    tz_cellobjfeats_clst(combclass,combcellidx,combobj,post,ncluster)
%TZ_CELLOBJFEATS_CLST Cell-level features for clustered objects
%   FEATURES = 
%       TZ_CELLOBJFEATS_CLST(COMBCLASS,COMBCELLIDX,COMBOBJ,POST,NCLUSTER)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function features=tz_cellobjfeats_clst(combclass,combcellidx,combobj,post,ncluster)
%OVERVIEW:
%   calculate cell-level features from clustered objects
%PARAMETERS:
%   combclass - class index
%   combcellidx - cell index
%   combobj - object-level features
%   post - clusters
%   ncluster - number of clusters
%RETURN:
%   features - cell-level features MCF
%DESCRIPTION:
%
%HISTORY:
%   27-MAY-2004 Initial write TINGZ
%

if size(post,2)>1
    [m,objidcs]=max(post');
    objidcs=objidcs';
else
    objidcs=post;
end

caclass=tz_findclass(combclass);

nclass=length(caclass);

if ~exist('ncluster','var') | ncluster==-1
   ncluster=max(objidcs);
end

% for i=1:nclass
%     ncells(i)=max(combcellidx(combclass==i));
% end

features={};
for i=1:ncluster
    subfeatures=tz_cellobjfeats(combclass(objidcs==i),combcellidx(objidcs==i),combobj(objidcs==i,:),nclass,[],{});
    features=tz_combinemcf(features,subfeatures);
end