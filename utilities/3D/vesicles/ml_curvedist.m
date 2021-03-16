function dists = ml_curvedist(pts,curve)
%ML_CURVEDIST Distance from a point to a curve.
%   DISTS = ML_CURVEDIST(PTS,CURVE) returns a vector of distances between
%   the [point array] PTS and [curve] CURVE. If the [curve] is closed, the
%   distance is postive if the point is outside the enclosed region of the
%   curve, othewise it is negative. In this function both PTS and CURVES
%   should be integer values.
%   
%   See also

%   20-Jan-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

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

if nargin < 2
    error('Exactly 2 arguments are required')
end

box = ml_boundbox([pts;curve]);
imageSize = box(2,:);

curvePoints = ml_showpts_2d(curve,'ln',0);
curveImage = ml_obj2img(curvePoints,imageSize);
distanceMap = bwdist(curveImage);

%If the curve is closed, we assign signs.
if all(curve(1,:)==curve(end,:))
    curveImage = imfill(curveImage,'hole');
    curveImage(curveImage==1) = -1;
    curveImage(curveImage==0) = 1;
    distanceMap = distanceMap.*curveImage;
end

pointIndices = sub2ind(imageSize,pts(:,1),pts(:,2));

dists = distanceMap(pointIndices);


