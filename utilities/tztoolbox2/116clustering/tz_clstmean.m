function y = tz_clstmean(data,clstlabel)
%TZ_CLSTMEAN Mean of features in clusters.
%   Y = TZ_CLSTMEAN(DATA,CLSTLABEL) returns the mean of data in each
%   clusters in CLSTLABEL, which is the vector of cluster labels of DATA.
%   The number of clusters is the maximum number in CLSTLABEL. DATA is a
%   Nxp matrix for N samples and P features. So Y is a PxK row vector.
%   The first P elements in Y is the mean of data in cluster 1, and so 
%   on. A better way to do this is to use the function:
%   Y = TZ_CLSTFEAT(DATA,CLSTLABEL,{'mean'});
%
%   See also TZ_CLSTFEAT

%   22-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

if size(data,1)~=length(clstlabel)
    error('the size of data should be same as the number of labels.');
end

ncluster = max(clstlabel);
nfeature = size(data,2);
y=[];

for i=1:ncluster
    clusterIdx = find(clstlabel==i);
    if isempty(clusterIdx)
        y = [y zeros(1,nfeature)];
    else
        y = [y mean(data(clusterIdx,:),1)];
    end
    
end
