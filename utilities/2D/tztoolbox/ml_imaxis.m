function [imgaxis,axln,dists,borders]=ml_imaxis(img)
%ML_IMAXIS Extract and process medial axis from an image.
%   IMGAXIS = ML_IMAXIS(IMG) returns an image is superimposition of IMG and
%   its mecial axis. This function will take all pixels above intensity 0 
%   as objects.
%   
%   [IMGAXIS,AXLN,DISTS,BORDERS] = ML_IMAXIS(...) also returns extracted
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


% Copyright (C) 2007  Murphy Lab
% Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.
%
% For additional information visit http://murphylab.web.cmu.edu or
% send email to murphy@cmu.edu

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
