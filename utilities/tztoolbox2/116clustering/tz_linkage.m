function Z = tz_linkage(X,distfun,groupmeth) 
%TZ_LINKAGE Hierarchical clustering information.
%   Z = TZ_LINKAGE(X,DISTFUN,GROUPMETH) computes the hierachical
%   clustering information on the data matrix X. DISTFUN is a string
%   specifying distance function:
%      'euclid'    --- Euclidean metric
%      'seuclid'   --- Standardized Euclid metric
%      'cityblock' --- City Block metric
%      'mahal'     --- Mahalanobis metric
%      'minkowski' --- Minkowski metric
%      'pvalue'    --- p-value distance
%   GROUPMETH is a string specifying group method:
%      'avgroup' 
%      'single'
%      'complete'
%      'average'
%      'centroid' 
%      'ward'
%
%   See also

%   ??-???-???? Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 3
    error('Exactly 3 arguments are required')
end

switch groupmeth
    case 'avgroup'
        
        Q=X;
        clst={};
        for i=1:size(X,1)
            clst{i}=X(i,:);
        end
        Z=[];
        
        for i=1:(size(X,1)-1)
            Y=tz_pdist(Q,distfun);
            [mindist,index]=min(Y);
            [m,n]=size(Y);
            m = (1 + sqrt(1+8*n))/2;
            
            S = zeros(m);
            I = ones(n,1);
            J = [1 (m-1):-1:2]';
            I(cumsum(J)) = (2:m);
            I2=cumsum(I);
            S(I2(index)) = 1;
            [posi,posj]=find(S==1);
            clst{end+1}=[clst{posi};clst{posj}];
            Q(posi,:)=NaN;
            Q(posj,:)=NaN;
            Q=[Q;mean(clst{end})];
            Z(i,:)=[posi posj mindist*10];
        end
    otherwise
        Y=tz_pdist(X,distfun);
        Z=linkage(Y,groupmeth);
end
