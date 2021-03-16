function tz_plotgraph(V,E,plotdim,pointstyle,linestyle,plotwndsize)
%TZ_PLOTGRAPH Plot a 2D undirected graph.
%   TZ_PLOTGRAPH(V,E,PLOTDIM,POINTSTYLE,LINESTYLE,PLOTWNDSIZE) plot a 2D
%   undirected graph for the graph (V,E), where V is a N-colunm matrix of
%   N-D nodes and E is a 2-column matrix of edges. PLOTDIM is a vector
%   specifying the dimensions for plotting. It can only be length 2 or 3. 
%   POINT specifies point plotting style in PLOT function, like 'ro' for
%   red circles.
%   LINESTYLE is the parameter for 'linestyle' option in PLOT function.
%   PLOTWNDSIZE is a 2-row vector giving boundray of the plotting window.

%   28-MAY-2004 Initial write T. Zhao
%   03-JUN-2004 Modified T. Zhao
%       -add plotwndsize
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 5
    error('Exactly 5 arguments are required')
end

if isempty(plotdim)
    plotdim=1:2;
end

if length(plotdim)~=2 & length(plotdim)~=3
    error('dimension error');
end

nv=size(V,1);
ne=size(E,1);

if(length(plotdim)==2)
    if ~isempty(pointstyle)
        plot(V(:,plotdim(1)),V(:,plotdim(2)),pointstyle);
    end
    if ~isempty(plotwndsize)
        if(size(V,1)>1)
            range=[min(V(:,plotdim));max(V(:,plotdim))];
            center=sum(range,1)/2;
            plotwnd=[center-plotwndsize/2;center+plotwndsize/2];
            xlim(plotwnd(:,1));
            ylim(plotwnd(:,2));
        end
    end
else
    plot3(V(:,plotdim(1)),V(:,plotdim(2)),V(:,plotdim(3)),pointstyle);
    if ~isempty(plotwndsize)
        if(size(V,1)>1)
            range=[min(V(:,plotdim));max(V(:,plotdim))];
            center=sum(range,1)/2;
            plotwnd=[center-plotwndsize/2;center+plotwndsize/2];
            xlim(plotwnd(:,1));
            ylim(plotwnd(:,2));
            zlim(plotwnd(:,3));
        end
    end
end

for i=1:ne
    if(length(plotdim)==2)
        line([V(E(i,1:2),plotdim(1))],[V(E(i,1:2),plotdim(2))], ...
            'Color',[1,0,0],'LineStyle',linestyle);
    else
        line([V(E(i,1:2),plotdim(1))],[V(E(i,1:2),plotdim(2))], ...
            [V(E(i,1:2),plotdim(3))],'Color',[0.8,0,0], ...
            'LineStyle',linestyle);
    end
end

