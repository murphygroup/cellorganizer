function a2=tz_normangle(a,option)
%TZ_NORMANGLE Obsolete. See ML_NORMANGLE.
%   A2 = TZ_NORMANGLE(A,OPTION) returns the normalized angle of A according 
%   to OPTION:
%       '180' - [0 180)
%       '360' - [0,360)
%       '90 - [0,90]
%       'a2r' - angle to radian
%       'r2a' - radian to angle

%   ??-OCT-2004 Initial write T. Zhao
%   01-NOV-2004 Modified T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

error(tz_genmsg('of','tz_normangle','ml_normangle'));

if nargin < 2
    error('Exactly 2 arguments are required')
end

switch option
    case '180'  %[0 180)
        a2 = tz_normangle(a,'360');
        a2(a2>180) = 360-a2(a2>180);
    case '360' 
        a2=mod(a,360);  %[0,360)
    case '90'   %[0,90]
        a2=abs(a);
        a2=mod(a2,180);
        ri=find(a2>=90);
        if ~isempty(ri)
            a2(ri)=180-a2(ri);
        end
    case 'a2r'  %angle to radian
        a2=a*pi/180;
    case 'r2a'  %radian to angle
        a2=a*180/pi;
end
