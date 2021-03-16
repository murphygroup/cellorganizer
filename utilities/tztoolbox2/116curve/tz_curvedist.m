function dists = tz_curvedist(pts,curve)
%TZ_CURVEDIST Obsolete. See ML_CURVEDIST.
%   DISTS = TZ_CURVEDIST(PTS,CURVE) returns a vector of distances between
%   the [point array] PTS and [curve] CURVE. If the [curve] is closed, the
%   distance is postive if the point is outside the enclosed region of the
%   curve, othewise it is negative. In this function both PTS and CURVES
%   should be integer values.
%   
%   See also

%   20-Jan-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

error(tz_genmsg('tz_curvedist','ml_curvedist'));

if nargin < 2
    error('Exactly 2 arguments are required')
end

box = tz_boundbox([pts;curve]);
imageSize = box(2,:);

curvePoints = tz_showpts_2d(curve,'ln',0);
curveImage = tz_obj2img(curvePoints,imageSize);
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


