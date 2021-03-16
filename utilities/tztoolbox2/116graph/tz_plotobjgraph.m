function tz_plotobjgraph(objcoords,procdim,plotdim,option,plotwndsize)
%TZ_PLOTOBJGRAPH Plot object graph.
%   TZ_PLOTOBJGRAPH(OBJCOORDS,PROCDIM,PLOTDIM,OPTION) plots
%   objects as a graph. OBJCOORDS is a matrix of coordinates of object
%   centers. If OPTION is 'null', the graph has no edges. If OPTIONS is
%   'mst', the graph will be a minimal spanning tree (MST) based on
%   dimensions specified by PROCDIM, which is a vector of column
%   indices of OBJCOORDS. Each row of OBJCOORDS represents an object.
%   
%   TZ_PLOTOBJGRAPH(OBJCOORDS,PROCDIM,PLOTDIM,OPTION,PLOTWNDSIZE)
%   also specifies the size of the plotting window.
%
%   Both PLOTDIM and PLOTWNDSIZE are parameters for TZ_PLOTGRAPH. See
%   TZ_PLOTGRAPH for more details

%   27-MAY-2004 Initial write T. Zhao
%   03-JUN-2004 Modified T. Zhao
%       -add plotwndsize
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 3
    error('3 or 4 arguments are required')
end

if ~exist('plotwndsize','var')
    plotwndsize=[];
end

if isempty(plotdim)
    plotdim=1:2;
end

if isempty(procdim)
    V=objcoords;
else
    V=objcoords(:,procdim);
end

switch option
case 'null'
    E=[];
case 'mst'
    if(size(V,1)<=1)
        E=[];
    else
        E=tz_mstree(V,'eu',1);
        E=E';
    end
otherwise
    warning('invalid option');
end

tz_plotgraph(V,E,plotdim,'x','-',plotwndsize);