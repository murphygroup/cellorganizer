function y = tz_clstfeat(data,clstlabel,featset)
%TZ_CLSTFEAT Summarized features on clusters.
%   Y = TZ_CLSTFEAT(DATA,CLSTLABEL,FEATSET) returns a vector of summarized
%   features of DATA on clusters, which are repesented by the vector 
%   CLSTLABEL. The number of clusters is the maximum number in CLSTLABEL.
%   DATA is a NxP matrix for N samples and P features. FEATSET is a cell
%   array of strings, which specify hwo to summarize DATA. A string can be
%       'mean' - mean value
%       'var' - variance
%       'std' - standard deviance
%   So Y is a row vector and it has PxMxK rows if there are K clusters and 
%   FEATSET contains M strings.
%   The first P elements in Y contain the first summarized set of data in 
%   cluster 1. The second P elements contain the second summarized set of
%   data in cluster 1 if there is a second summarized set. If not, it is 
%   the the first summarized set of data in cluster 2, and so on.
%   It is always like:
%   [cluster1(featset1 featset2 ...) cluster2(featset1 featset2 ...) ...]

%   22-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 3
    error('Exactly 3 arguments are required')
end

if size(data,1)~=length(clstlabel)
    error('the size of data should be same as the number of labels.');
end

ncluster = max(clstlabel);
nfeature = size(data,2);
nfeatset = length(featset);
y=[];

for i=1:ncluster
    clusterIdx = find(clstlabel==i);
    subdata = data(clusterIdx,:);
    for j=1:length(featset)
        if isempty(clusterIdx)
            y = [y zeros(1,nfeature)];
        else
            switch(featset{j})
            case 'mean'
                y = [y mean(subdata,1)];
            case 'var'
                if(length(clusterIdx)==1)
                    y = [y zeros(1,nfeature)];
                else
                    y = [y var(subdata)];
                end
            case 'std'
                if(length(clusterIdx)==1)
                    y = [y zeros(1,nfeature)];
                else
                    y = [y std(subdata)];
                end        
            end
        end
    end
end
