function tz_buildtree(X,labels,distfun,groupmeth,iszscore)
%TZ_BUILDTREE Under construction.

% function tz_buildtree(X,labels,distfun,groupmeth,iszscore)
%OVERVIEW:
%   Build a hierachical tree for data X
%PARAMETERS:
%   X - data for clustering
%   labels - names for data
%   distfun - distance function
%   groupmeth - grouping method
%   iszscore - zscore or not
%RETURN:
%   no return
%DESCRIPTION:
%   This function shows a dendrogram for the data without any return value
%
%HISTORY:
%   ??-???-???? Initial write TINGZ

if iszscore==1
    X=zscore(X);
end

Z=tz_linkage(X,distfun,groupmeth);
xc_dendrogram(Z,0,labels);
