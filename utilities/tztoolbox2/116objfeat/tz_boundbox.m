function [box,s] = tz_boundbox(pts)
%TZ_BOUNDBOX Obsolete. See ML_BOUNDBOX.
%   BOX = TZ_BOUNDBOX(PTS) returns a 2x2 matrix which is the bounding box 
%   of the [point array] PTS. The first row of BOX is the left top corner
%   coordinate and the second row of BOX is the right bottom corner
%   coordinate. Here left,top,right,bottom is based on [image coordinate
%   system].
%
%   [BOX,S] = TZ_BOUNDBOX(PTS) also returns the size of the bounding box.
%   S(1) is the height and S(2) is the width.
%
%   See also

%   20-Jan-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

error(tz_genmsg('of','tz_boundbox','ml_boundbox'));

if nargin < 1
    error('Exactly 1 argument is required')
end

topleft = min(pts,[],1);
bottomright = max(pts,[],1);

box = [topleft; bottomright];
s = box(2,:)-box(1,:)+[1 1];

