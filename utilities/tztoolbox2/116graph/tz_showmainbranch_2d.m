function tz_showmainbranch_2d(objcoords,mst,ptstyle,lnstyle)
%TZ_SHOWMAINBRANCH_2D Show the main branch of MST.
%   TZ_SHOWMAINBRANCH_2D(OBJCOORDS,MST,POINTSTYLE,LNSTYLE) shows the main
%   branch of points of the minimal spanning tree of nodes OBJCOORDS.
%   POINTSTYPE and LINESTYLE are parameters for TZ_PLOTGRAPH.. MST is the
%   calculated minimal spanning tree. If MST is empty, the tree will be
%   extracted from OBJCOORDS automatically.

%   ??-???-???? Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University


if nargin < 4
    error('Exactly 4 arguments are required')
end

if isempty(mst)
    mst=tz_mstree(objcoords(:,1:2),'eu',1);
end

maxpath=tz_treemainbranch(mst');

tz_plotptser_2d(objcoords(maxpath,1:2),'x','-')
