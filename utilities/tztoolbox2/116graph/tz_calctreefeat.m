function [feature,name,degree]=tz_calctreefeat(tree,V)
%TZ_CALCTREEFEAT Calcuate features of a tree.
%   FEATURE = TZ_CALCTREEFEAT(TREE) returns a row vector containing
%   features of the tree TREE, which is an Nx2 matrix if there are N
%   edges. Each element in TREE is the label of a node.
%   
%   FEATURE = TZ_CALCTREEFEAT(TREE,V) works on a tree with 2D coordinated
%   nodes V, which has coordinates for a node in each row.
%   
%   [FEATURE,NAME,DEGREE] = TZ_CALCTREEFEAT(...) also returns the names
%   of the features and degrees of all the nodes.
%
%   The features are:
%       'mean of branch length', 'variance of branch length', 
%       'longest branch', 'shortest branch', 'max degree', 'leave numbers',
%       'degree mean', 'degree var', 'length of longest trunk',
%       'node number of longest trunk', 'mean of angle', 
%       'variance of angle', 'mean of degree of main branch',
%        'variance of degree of main branch'

%   ??-???-2004 Initial write  Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('1 or 2 arguments are required')
end

if ~exist('V','var')
    V=[];
end

name={'mean of branch length','variance of branch length', ...
        'longest branch','shortest branch','max degree'...
        'leave numbers','degree mean','degree var', ...
        'length of longest trunk','node number of longest trunk',...
        'mean of angle','variance of angle', ...
        'mean of degree of main branch', ...
        'variance of degree of main branch'};

if isempty(tree)
    feature=zeros(1,14);
    return;
end

feature(1)=mean(tree(:,3));
feature(2)=var(tree(:,3));
feature(3)=max(tree(:,3));
feature(4)=min(tree(:,3));

nnode=max([tree(:,1);tree(:,2)]);

for i=1:nnode
    degree(i)=sum(tree(:,1)==i)+sum(tree(:,2)==i);
end

feature(5)=max(degree);
feature(6)=sum(degree==1);
feature(7)=mean(degree);
feature(8)=var(degree);
[path,feature(9)]=tz_treemainbranch2(tree);
feature(10)=length(path);

if isempty(V)
    feature(11)=0;
    feature(12)=0;
else
    a=0;
    for i=1:feature(10)-2
        a(i)=tz_findangle_2d(V(path(i),1:2),V(path(i+1),1:2),V(path(i+2),1:2));
    end
    feature(11)=mean(a);
    feature(12)=var(a);
end

for i=1:feature(10)
    mdeg(i)=degree(path(i));
end

feature(13)=mean(mdeg);
feature(14)=var(mdeg);

