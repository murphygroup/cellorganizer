function [path,len]=tz_treepath(snode,tnode,tree,option)
%TZ_TREEPATH Find a path in a tree.
%   PATH = TZ_TREEPATH(SNODE,TNODE,TREE) return the path from the node
%   SNODE to TNODE in the tree TREE. Each node in a tree is represented 
%   by a unique natrual number. PATH is a vector of nodes.
%   
%   PATH = TZ_TREEPATH(SNODE,TNODE,TREE,OPTION) specifies the form of
%   TREE by OPTION. If option is 1, TREE is has matrix form with wights
%   at the corresponding rows and columns. Otherwise, TREE has a
%   [node, node, weight] form.
%   
%   [PATH,LEN] = TZ_TREEPATH(...) also return the length of the path.

%   ??-???-???? Initial write T. Zhao
%   05-NOV-2004 Modified T. Zhao
%       - add comments
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 3
    error('3 or 4 arguments are required')
end

if ~exist('option','var')
    option=1;
end

if(option==1)
    m=tz_tree2matrix(tree);
    
    m=tril(m)+triu(m)';
    m=m+m';
else
    m=tree;
end


path=[snode,snode];
len=0;

curnode=snode;

while(curnode~=tnode)
    cnode=find(m(curnode,:)~=0);
    if any(cnode==path(end-1))
        cnode(cnode==path(end-1))=[];
    end
    
    if isempty(cnode)
        path(end)=[];
        m(curnode,:)=0;
        m(:,curnode)=0;
        curnode=path(end);
    else
        curnode=cnode(1);
        path=[path,curnode];
    end
end

path(1)=[];

for i=1:length(path)-1
    len=len+m(path(i),path(i+1));
end

    