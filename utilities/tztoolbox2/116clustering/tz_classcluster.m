function labels=tz_classcluster(combobj,combclass,k,norm,t)
%TZ_CLASSCLUSTER Classwise clustering.
%   LABEL = TZ_CLASSCLUSTER(COMBOBJ,COMBCLASS,K,NORM) clusters all objects
%   based on their feature matrix COMBOBJ class by class. COMBCLASS is a
%   vector of class labels. COMBOBJ and COMBCLASS have the same number
%   of rows. K is the cluster number of each class. If K is -1, the number
%   of clusters will be determined automatically. But it only searches K
%   from 1 to 15. NORM is the way of normalizing data:
%       0 - no normalization
%       1 - z-score upon all data
%       2 - z-score class by class
%   
%   LABEL = TZ_CLASSCLUSTER(COMBOBJ,COMBCLASS,K,NORM,T) also lets users
%   specify the way of searching the optimal number of clusters. T is 
%   a vector. The function will try different number of K in T(1:end-1)
%   for T(end) times and find the best one.

%   ??-???-???? Initial write T. Zhao
%   19-MAY-2004 Modified T. Zhao
%       - add parameter t
%   30-Aug-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 4
    error('4 or 5 arguments are required')
end

if ~exist('t','var')
    t = [];
end
if isempty(t)
    t=[1:30 10];
end

nclass=max(combclass);

fnorm=combobj;
if norm==1
    fnorm=zscore(fnorm);
end

labels=zeros(size(combclass));
p=size(combobj,2);
allcenters=[];
for i=1:nclass
    i
    subfnorm=fnorm(combclass==i,:);
    if norm==2
        subfnorm=zscore(subfnorm);
    end
    
    if(k>0)
        randidx = randperm(size(subfnorm,1));
        seeds = subfnorm(randidx(1:k),:);
        
        % Do k-means clustering
        options = zeros(1,14);
%         options(14)=100;
        
        [centers,options2,post,errlog]=tz_kmeans(seeds,subfnorm,options,-1);
    else
        if k==0
            [aics,centers,posts]=tz_kmeansaic(subfnorm,t(1:end-1),t(end),0,0);
        else
            [aics,centers,posts]=tz_kmeansaic(subfnorm,t(1:end-1),t(end),0,-1);
        end
        [c,m,n]=ml_min(aics);
        post=posts{m,n};
        allcenters=[allcenters;centers{m,n}];
    end
    
    [tmp,post]=max(post');
    ncluster=max(post);
    
    labels(combclass==i)=post'+max(labels);
end

if k==-2
    % Do k-means clustering
    options = zeros(1,14);
    [newcenters,options2,post,errlog]=tz_kmeans(allcenters,fnorm,options,-1);
    [tmp,labels]=max(post');
    labels=labels';
end
