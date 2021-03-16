function as = tz_ptlnangle(pts,lnpts)
%TZ_PTLNANGLE Calculate angles between points and a line.
%   AS = TZ_PTLNANGLE(PTS,LNPTS) returns the angles between the array of
%   points PTS and the line LNPTS. LNPTS could be 2x2 matrix or 1x3 vector,
%   corresponding to different form of lines. 2x2 matrix contains
%   [start_point;end_point] and 1x3 vector contains [start_point, angle].
%   
%   See also TZ_PTLNDIST

%   24-Apr-2005 Initial write T T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

x0=pts(:,1);
y0=pts(:,2);
x1=lnpts(1,1);
y1=lnpts(1,2);
if length(lnpts)==3
    ra=lnpts(3)*pi/180;
    x2=x1+cos(ra);
    y2=y1+sin(ra);
else
    x2=lnpts(2,1);
    y2=lnpts(2,2);
end

v1=[x2-x1,y2-y1];
as=acos(pts*v1'./sqrt(sum(pts.^2,2))/sqrt(sum(v1.^2)));
as=as.*sign(pts*[-v1(2); v1(1)]);