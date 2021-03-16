function [post2,rcidx,rridx]=ml_stdclst(post,threshold)
%ML_STDCLST Remove small clusters.
%   POST2 = ML_STDCLST(POST) removes empty clusters in POST and returns
%   the new clusters POST2. Both POST and POST2 are [label matrix].
%   
%   POST2 = ML_STDCLST(POST,THRESHOLD) removes clusters with size no 
%   greater than THRESHOLD in POST. Samples not belonging to any left
%   clusters will also be removed.
%
%   [POST2,RCIDX,RRIDX] = ML_STDCLST(...) also returns the indices of 
%   removed clusters RCIDX and removed samples RRIDX.
%
%   Example:
%       if x = [1 0 0 0    then ml_stdclst(post) is [1 0 0
%               1 0 0 0                              1 0 0
%               0 1 0 0                              0 1 0
%               0 0 0 1]                             0 0 1]
%
%       and [post2,rcidx,rridx]=ml_stdclst(post,1) will let post2 be
%       [1 1]', rcidx be [2 3 4] and rridx be [3 4]'.

%   ??-???-???? Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('1 or 2 arguments are required')
end

if nargin < 2
    threshold = 0;
end

if ~all( ismember(post(:),[0,1]) )
    error('The first argument must contain only 0 or 1');
end

if any(sum(post,2) ~= 1)
    error('Invalid coding. Sum of columns must be one');
end

clstsize=sum(post,1);

post2=post;
rcidx=find(clstsize<=threshold);
post2(:,rcidx)=[];

rridx=find(sum(post2,2)==0);
post2(rridx,:)=[];
