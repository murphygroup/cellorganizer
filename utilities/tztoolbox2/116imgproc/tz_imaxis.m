function [imgaxis,axln,dists,borders]=tz_imaxis(img)
%TZ_IMAXIS Obsolete. See ML_IMAXIS.
%   IMGAXIS = TZ_IMAXIS(IMG) returns an image is superimposition of IMG and
%   its mecial axis. This function will take all pixels above intensity 0 
%   as objects.
%   
%   [IMGAXIS,AXLN,DISTS,BORDERS] = TZ_IMAXIS(...) also returns extracted
%   data from the medial axis represenation. AXLN is the coordinates of
%   medial axis, in which each row has X and Y coordinates. DISTS is the
%   width along medial axis. BORDERS is the coordinates of the contour of
%   the medial axis representation. PTS is the coordinate representation of
%   the shape.
%   
%   See also

%   ??-???-2004 Initial write T. Zhao
%   04-NOV-2004 Modified TINGZ
%       - add comments
%   Copyright (c) Murphy Lab, Carnegie Mellon University

error(tz_genmsg('of','tz_imaxis','ml_axis'));

[edgex,edgey]=find(img>0);

imgaxis=img;
minx=min(edgex(:));
maxx=max(edgex(:));
k=1;
for i=minx:maxx
    curedgey=edgey(edgex==i);
    if ~isempty(curedgey)
        maxy=max(curedgey(:));
        miny=min(curedgey(:));
        dists(k)=maxy-miny+1;
        borders(k,:)=[miny,maxy];
        axln(k,:)=[i,round((maxy+miny)/2)];
        imgaxis(axln(k,1),axln(k,2))=2;
        k=k+1;
    end
end
