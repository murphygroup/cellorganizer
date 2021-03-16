function a=tz_findangle_2d(pt1,pt2,pt3)
%TZ_FINDANGLE_2D Calculate angles according to three points.
%   A = TZ_FINDANGLE_2D(PT1,PT2,PT3) returns the angle of two lines, which
%   have ends PT1,PT2 and ends PT2,PT3. PT1, PT2 and PT3 are row vectors
%   with two elements. A has unit radian and range [0,2*pi). If the line 
%   with ends PT2,PT1 rotates A degree counterwisely around PT2, it will be
%   overlapped with the line with ends PT2,PT3.
%   
%   See also

%   16-Sep-2005 Initial write T. Zhao
%   09-Feb-2005 Modified T. Zhao
%       - Fix 0 angle bug
%   Copyright (c) Murphy Lab, Carnegie Mellon University

v1=pt1-pt2;
v2=pt3-pt2;

if all(v1==0) | all(v2==0)
    a=0;
    warning('0 vectors');
    return;
end

a=acos(sum(v1.*v2)/sqrt(sum(v1.*v1))/sqrt(sum(v2.*v2)));

switch sign(sum(v2.*[-v1(2) v1(1)]))
case 1
    a=a;
case -1
    a=2*pi-a;
end

