function ds = tz_ptlndist(pts,lnpts)
%TZ_PTLNDIST Calculate distances between points and a line.
%   DS = TZ_PTLNDIST(PTS,LNPTS) returns the distances between the points
%   PTS and the line LNPTS. LNPTS could be 2x2 matrix or 1x3 vector,
%   corresponding to different form of lines. 2x2 matrix contains
%   [start point;end points and 1x3 vector contains [start point, angle].
%   
%   See also TZ_PTLNANGLE

%   24-Apr-2005 Initial write  T. Zhao
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
ds=abs((x2-x1)*(y1-y0)-(x1-x0)*(y2-y1))/sqrt((x2-x1)^2+(y2-y1)^2);

