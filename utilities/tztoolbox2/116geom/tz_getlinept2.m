function pts=tz_getlinept2(s,a,len,isend)
%TZ_GETLINEPT2 Obsolsete. See ML_GETLINEPT2.
%   TZ_GETLINEPT2(S,A,LEN) returns coordinates of points on a line segment
%   with starting point S, angle A and length LEN.
%   
%   TZ_GETLINEPT2(S,A,LEN,ISEND) only returns the two ends of the line
%   segment if ISEND is 1. Otherwise, it is the same as 
%   TZ_GETLINEPT2(S,A,LEN)

%   ??-???-???? Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

error(tz_genmsg('of','tz_getlinept2','ml_getlinept2'));

if ~exist('isend','var')
    isend=0;
end

a=mod(a,360);
len=round(len);

switch a
case 0
    t=[s(1)+len,s(2)];
case 90
    t=[s(1),s(2)+len];
case 180
    t=[s(1)-len,s(2)];
case 270
    t=[s(1),s(2)-len];
otherwise
    ra=a*pi/180;
    t=round([s(1)+cos(ra)*len,s(2)+sin(ra)*len]);
end

if isend
    pts=[s;t];
else
    pts=tz_getlinept(s,t);
end