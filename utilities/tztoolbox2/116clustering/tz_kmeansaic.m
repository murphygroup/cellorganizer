function [aics,centers,posts]=tz_kmeansaic(combobj,ks,tr,norm, ...
    rmsmall,rtlabel)
%TZ_KMEANSAIC K-means clustering with AIC.
%   AICS = TZ_KMEANSAIC(X,KS,TR,NORM,RMSMALL) returns AIC values
%   of the clusters of the data in the matrix X. The clusters were found by
%   batch k-means on the rows of the data. AICS is matrix in which each 
%   row vector contains AIC values from different numbers of clusters
%   KS and each column vector contains AIC values from different trails
%   of TR times. If NORM is 0, the original data will be used, or X
%   will be z-scored before clustering. RMSMALL determines the allowable 
%   size of the smallest cluster. The sizes of final clusters will be all
%   greater than RMSMALL.
%   
%   [AICS,CENTERS,POSTS] = TZ_KMEANSAIC(...) will also returns the 
%   centers of clusters and group membership for each sample. Both
%   CENTERS and POSTS are cell array, with the same dimensions as
%   those of AICS.
%   
%   [AICS,CENTERS,POSTS] = 
%       TZ_KMEANSAIC(X,KS,TR,NORM,RMSMALL,RTLABEL) lets users
%   specify the format of POSTS. RTLABEL =1 returns scalar labels,
%   otherwise POSTS is a martix with 1-N coding.

%   17-MAY-2004 Initial write T. Zhao
%   06-SEP-2004 Modified T. Zhao
%       - Add the parameter rmsmall
%   30-Aug-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

warning('Seeding schema changed');

if ~exist('rtlabel','var')
    rtlabel=0;
end

i=1;

%data size
nobj=size(combobj,1);

%0: no preprocessing
if norm==0
    fnorm=combobj;
else
    fnorm=zscore(combobj);
end
    

%trials loop
for t=1:tr
    [num2str(t) 'trial']
    i=1;
%     randidx=randperm(nobj);
    options = zeros(1,14);
    %cluster numbers loop
    for k=ks
        if k==1
            post=ones(nobj,1);
            centers{i,t}=mean(fnorm,1);
        else
            seeds=num2str(k);    
%             if rmsmall==1
            [centers{i,t},options,post,errlog] = ...
                tz_kmeans(seeds,fnorm,options,rmsmall);
%             else
%                 [centers{i,t},options,post,errlog]=kmeans(seeds,fnorm,options);
%             end
        end
        [aics(i,t), clusterinfo] = ...
            tz_aicbic(fnorm,ml_post2label(post),'euc');
        
        if rtlabel==1
            posts{i,t}=ml_post2label(post);
        else
            posts{i,t}=post;
        end
        
        i=i+1;
    end
end