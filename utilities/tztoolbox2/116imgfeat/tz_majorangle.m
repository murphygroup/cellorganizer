function theta=tz_majorangle(img)
%TZ_MAJORANGLE Obsolete. See ML_MAJORANGLE

%UUU
%function theta=tz_majorangle(img)
%
%OVERVIEW:
%   find angle of the major axis of an image
%PARAMETERS:
%   img - input image
%RETURN:
%   theta - angle
%DESCRIPTION:
%
%HISTORY:
%   ??-???-2004 Initial write TINGZ
%   01-NOV-2004 Modified TINGZ
%   24-Mar-2005 Modified TINGZ
%       - debugged

error(tz_genmsg('of','tz_majorangle','ml_majorangle'));

mom = tz_moment2(img);
weights=img(find(img>0));

center=[mom.cx,mom.cy];

theta = .5 * atan((mom.mu02 - mom.mu20)/2/mom.mu11)+sign(mom.mu11)*pi/4;

ntheta=[cos(theta),sin(theta)];
[x,y]=find(img>0);
x=x-mom.cx;
y=y-mom.cy;
imgskew=tz_wmoment([x,y]*ntheta',weights,3);

if imgskew<0
    theta=theta+pi;
end