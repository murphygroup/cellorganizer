function theta = tz_ptsmangle(pts)
%TZ_PTSMANGLE Calculate the orientation of point sets.
%   THETA = TZ_PTSMANGLE(PTS) returns the major angle of the array of
%   points PTS. The unit is radian.
%   
%   See also

%   24-Apr-2005 Initial write  T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argument is required')
end

center=mean(pts,1);
covmat=cov(pts(:,1),pts(:,2),1);
mu20=covmat(1,1);
mu02=covmat(2,2);
mu11=covmat(1,2);

theta = .5 * atan((mu02 - mu20)/2/mu11)+sign(mu11)*pi/4;

ntheta=[cos(theta),sin(theta)];

normpts=ml_addrow(pts,-center);
ds=tz_ptlndist(normpts,[0 0 theta*180/pi]);
imgskew=ml_wmoment(normpts*ntheta',ds,3);
% imgskew=mean(ds.*(normpts*ntheta').^3);
if imgskew<0
    theta=theta+pi;
end