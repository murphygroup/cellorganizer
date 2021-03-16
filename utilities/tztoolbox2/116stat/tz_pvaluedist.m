function dist = tz_pvaluedist(p1,p2,b)

%function dist = tz_pvaluedist(p1,p2,b)
%
%OVERVIEW:
%   calculate pvalue distance
%PARAMETERS:
%   p1 - vector 1
%   p2 - vector 2
%   b - calibration
%RETURN:
%   dist - calculated distance
%DESCRIPTION:
%   dist(p1,p2)=p1*p2-2*p1*p2
%
%HISTORY:
%   ??-???-???? Initial write TINGZ
%   24-DEC-2003 Modified TINGZ

if (~exist('b','var'))
    b=0.5;
end

%dist = tz_pvalue_dist(p1,p2)
ps1=p1;  %./(p1+b);
ps2=p2;  %./(p2+b);

dist = mean(ps1+ps2-2*ps1.*ps2);
