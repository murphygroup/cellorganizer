function pts2 = tz_moldshape(pts,dists)
%TZ_MOLDSHAPE Obsolete. See ML_MOLDSHAPE.
%   PTS2 = TZ_MOLDSHAPE(PTS,DISTS) returns an array of points that are form
%   the points PTS, which are relocated according to the distance
%   differences DISTS. PTS and DISTS must have the same number of rows.
%   
%   See also

%   11-Apr-2005 Initial write  T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

error(tz_genmsg('of','tz_moldshape','ml_moldshape'));

if nargin < 2
    error('Exactly 2 arguments are required')
end 

[tmp,s]=min(abs(pts(pts(:,1)>0,2)));

if s>1
    tpts=[pts(s:end,:);pts(1:s-1,:)];
else
    tpts=pts;
end
orgdists=sqrt(sum(tpts.^2,2));
pts2(:,1)=tpts(:,1).*(orgdists+dists)./orgdists;
pts2(:,2)=tpts(:,2).*(orgdists+dists)./orgdists;
