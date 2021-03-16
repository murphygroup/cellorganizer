function tree=tz_matrix2tree(m)
%TZ_MATRIX2TREE Convert a matrix to a graph
%   TREE = TZ_MATRIX2TREE(M) returns a graph, which is a three-column 
%   matrx. Each row of TREE represents an edge of the graph by the 
%   form [node 1, node 2, edge weight].
%   
%   See also TZ_TREE2MATRIX

%   ??-???-???? Initial write T. Zhao
%   04-NOV-2004 Modified T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argument is required')
end

tree=[];

if isempty(m)
    return;
end

nnode=size(m,1);

pnode = 1;

for i=1:nnode
    cnode=find(pnode(i,:)~=0);
    for j=1:length(cnode)
        tree=[tree;[i,cnode(j),m(i,cnode)]];
    end
end