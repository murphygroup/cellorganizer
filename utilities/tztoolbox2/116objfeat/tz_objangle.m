function ca=tz_objangle(objcofs,cof,ra)
%TZ_OBJANGLE Calculate angles of objects in a cell.
%   CA = TZ_OBJANGLE(OBJCOFS,COF,RA) returns a vector of angles of object
%   COFs OBJCOFS in a cell. OBJCOFS is a 2-column matrix and each row is 
%   a COF of an object. COF is a 1x2 vector representing the COF of the
%   nuclues. The cell has the orientation RA with the form 
%   [angle,fliplr,flipud]. The unit of angles is degree in [0,360).

%   ??-???-???? Initial write T. Zhao
%   30-OCT-2004 Modified T. Zhao
%       - add comments
%       - change function name: tz_cellangel-->tz_objangel
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 3
    error('Exactly 3 arguments are required')
end

nobj=size(objcofs,1);

objvec=[objcofs(:,1)-cof(1),objcofs(:,2)-cof(2)];
objangle=atan2(objvec(:,1),objvec(:,2));
ca=objangle*180/pi-ra(1);
if ra(2)==1
    ca=ca+180;
end

if ra(3)==1
    ca=-ca;
end

ca=tz_normangle(ca,'360');

