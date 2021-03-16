function m=tz_tree2matrix(tree)
%TZ_TREE2MATRIX Convert a weighted tree to a matrix.
%   M = TZ_TREE2MATRIX(TREE) returns a matrix representing a graph. The
%   value an element at the Ith row and Jth column in M is the weight of
%   edge of two nodes I and J. If there is no such an edge, the value will
%   be 0. Each row of TREE is like [node 1, node 2, edge weight].

%   ??-???-???? Initial write T. Zhao
%   05-NOV-2004 Modified T. Zhao
%       - add comments
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argument is required')
end

if isempty(tree)
    m=[];
    return;
end

nedge=size(tree,1);
m=zeros(nedge+1);

for i=1:nedge
    m(tree(i,1),tree(i,2))=tree(i,3);
end