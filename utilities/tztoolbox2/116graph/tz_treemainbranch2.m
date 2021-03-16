function [maxpath,maxlen]=tz_treemainbranch2(tree)
%TZ_TREEMAINBRANCH Find the main branch of a tree.
%   It is the same same as TZ_TREEMAINBRANCH, but it is faster.

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
tmpm=m;
cm=m~=0;

leaves=find(sum(cm,1)==1);
maxpath=[];
maxlen=0;

for i=1:length(leaves)-1
    snode=leaves(i);
    tnodes=leaves(i+1:end);
    
    path=[snode,snode];
    len=0;
    m=tmpm;
    curnode=snode;
    
    while(~isempty(tnodes))
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
        
        if any(curnode==tnodes)
            tnodes(tnodes==curnode)=[];
            len=0;
            for i=1:length(path)-1
                len=len+m(path(i),path(i+1));
            end
            
            if len>maxlen
                maxpath=path;
                maxlen=len;
            end
        end
    end
   
end
maxpath(1)=[];