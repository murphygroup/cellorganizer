function dists=tz_calcedgedist(V,E)
%TZ_CALCEDGEDIST Calcualte the lengths of all edges in a graph.
%   DISTS = TZ_CALCEDGEDIST(V,E) returns a vector of the lengths of all 
%   edges in the graph (V,E), where V is the set of nodes and E is the set
%   of edges.

%   ??-???-???? Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

dists=sqrt(sum((V(E(:,1),:)-V(E(:,2),:)).^2,2));
