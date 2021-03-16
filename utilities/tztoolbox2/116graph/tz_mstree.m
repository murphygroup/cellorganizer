function v=tz_mstree(X,s,t)
%TZ_MSTREE Minimal spanning tree of a weighted graph.
%   TREE = TZ_MSTREE(X,S,T) return a 3-ROW matrix representing the
%   minimal spanning tree (MST) of the data matrix X, which contains 
%   the positions of all points in the graph. Each row represents one point.
%   and there are N columns if it has N dimensions. S is a string 
%   specifying distance caluction method and T contains additional 
%   parameters for the method. See TZ_PDIST for more information. 
%   Each cloumn of TREE has the form [node 1;node 2;edge weights].
% 
%   See also TZ_MSTREE_L TZ_PDIST

%   04-NOV-2003 Initial write T. Zhao
%   05-NOV-2003 Modified T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 3
    error('Exactly 3 arguments are required')
end

if size(X,1)<=1
    v=[];
    return;
end

G=squareform(tz_pdist(X,s,t));

trilG=triu(NaN*ones(size(G)))+G;
symG=tril(G)+tril(G)';

pnum=size(G,1);
m=1:(pnum-1);
n=m;
c=m;

pset=1:pnum;
pnewset=[];

[c(1),m(1),n(1)]=tz_min(trilG);
pnewset=[pnewset,m(1),n(1)];
pset([m(1),n(1)])=[];
pnewset=sort(pnewset);
pset=sort(pset);

for i=2:(pnum-1)
    m1=symG(pnewset,:);
    m1(:,pnewset)=[];
    [c(i),mm,nn]=tz_min(m1);
    m(i)=pnewset(mm);
    n(i)=pset(nn);
    pnewset=[pnewset,n(i)];
    pset(nn)=[];   
    pnewset=sort(pnewset);
    pset=sort(pset);
end

v=[m;n;c];