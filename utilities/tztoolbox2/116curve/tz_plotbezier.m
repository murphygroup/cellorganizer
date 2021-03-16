function tz_plotbezier(ctrpts)
%TZ_BEZIER Plot 2D bezier curve.
%   Y = TZ_PLOTBEZIER(CTRPTS) plots 2D bezier defined by control points 
%   CTRPTS, which also determine the order of the bezier curve.
%
%   See also

%   27-Dec-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('Exactly 1 arguments are required')
end

t = 0:0.01:1;

x = tz_bernstein(t,ctrpts(:,1)');
y = tz_bernstein(t,ctrpts(:,2)');

plot(x,y);
hold on
% [ctrpts(:,1),idx] = sort(ctrpts(:,1));
% ctrpts(:,2) = ctrpts(idx,2);
plot(ctrpts(:,1),ctrpts(:,2),'ro');
plot(ctrpts(:,1),ctrpts(:,2),'r-.');
hold off


