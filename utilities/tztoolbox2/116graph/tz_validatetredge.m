function succ=tz_validatetredge(V,E,s,newv)
%TZ_VALIDATEEDGE Validate minimal spanning tree.
%   SUCC = TZ_VALIDATEEDGE(V,E,S,NEW) returns 1 if a node NEWV will not
%   change the structure of the minimal spanning tree E of nodes V. S is
%   a vector of nodes that are supposed to connected with NEWV. Otherwise,
%   it returns 0.

%   ??-???-???? Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

V=[V;newv];
newe=size(V,1);
newE=tz_mstree_l(V)';
ds=sqrt(sum((V(s,:)-newv).^2));
dists=tz_calcedgedist(V,[E;s size(V,1) 0])';
[E(:,3)' ds];
diff=abs(sum(newE(:,3))-sum(dists));
succ=diff<1e-3;