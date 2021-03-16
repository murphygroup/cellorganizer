function [maxpath,maxlen]=tz_treemainbranch(tree)
%TZ_TREEMAINBRANCH Find the main branch of a tree.
%   MAXPATH = TZ_TREEMAINBRANCH(TREE) returns the path of the main branch
%   of the tree TREE, which is a 3-column matrix. Each row of TREE is like 
%   [node 1, node 2, edge weight]. MAXPATH is a vector of natural numbers
%   representing nodes.
%   
%   [MAXPATH,MAXLEN] = TZ_TREEMAINBRANCH(...) also returns the length of
%   the path.

%   ??-???-???? Initial write T. Zhao
%   05-NOV-2004 Modified T. Zhao
%       - add comments
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argument is required')
end

m=tz_tree2matrix(tree);

m=tril(m)+triu(m)';
m=m+m';

cm=m~=0;

leaves=find(sum(cm,1)==1);

k=1;
for i=1:length(leaves)
    for j=i+1:length(leaves)
        if(leaves(i)>0 & leaves(j)>0)
            lp1=find(cm(leaves(i),:)==1);
            lp2=find(cm(leaves(j),:)==1);
            if(lp1==lp2)
                path{k}=[leaves(i),lp1,leaves(j)];
                len(k)=m(leaves(i),lp1)+m(leaves(j),lp1);
                if m(leaves(i),lp1)<=m(leaves(j),lp1)
                    leaves(i)=0;
                else
                    leaves(j)=0;
                end
            else            
                [path{k},len(k)]=tz_treepath(leaves(i),leaves(j),m,0);
            end
            k=k+1;
        end
    end
end

[maxlen,i]=max(len);
maxpath=path{i};