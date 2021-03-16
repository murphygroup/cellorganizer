function pts2 = tz_normapts(pts,center)
%TZ_NORMAPTS Let the points start with 0 angle.
%   PTS2 = TZ_NORMAPTS(PTS,CENTER) reorganize the array of points PTS so
%   that it starts from angle 0. Here the angle is the angle between a
%   boundary point and X axis.
%   
%   See also

%   11-Apr-2005 Initial write  T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

if isempty(center)
    center=mean(pts,1);
end
pts=ml_addrow(pts,-center);

%Exclude points at the negative side of X coordinate
tmppts=pts;
tmppts(tmppts(:,1)<0,2)=Inf;

%Find minimal Y value (THIS HAS A FLAW!)
[tmp,s]=min(abs(tmppts(:,2)));

if s>1
    pts2=[pts(s:end,:);pts(1:s-1,:)];
else
    pts2=pts;
end