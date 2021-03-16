function [a,t]=tz_alignshape(sh2,sh1,w)
%TZ_ALIGNSHAPE Obsolete. See ML_ALIGNSHAPE.
%   A = ML_ALIGNSHAPE(SH2,SH1) returns the rotation and scaling parameters
%   of transformation from shape SH2 to SH1. SH2 and SH1 should have the 
%   same size.
%   
%   A = ML_ALIGNSHAPE(SH2,SH1,W) takes account of weights for each point.
%
%   [A,T] = TZ_ALIGNSHAPE(...) also returns the translation vector.
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

error(tz_genmsg('of','tz_alignshape','ml_alignshape'));

if nargin < 2
    error('2 or 3 arguments are required')
end

if nargin<3
    w=[];
end

if isempty(w)
    w=ones(size(sh1,1),1);
end

X1=sum(w.*sh1(:,1));
X2=sum(w.*sh2(:,1));
Y1=sum(w.*sh1(:,2));
Y2=sum(w.*sh2(:,2));
W=sum(w);
Z=sum(w.*(sh2(:,1).^2+sh2(:,2).^2));
C1=sum(w.*sum(sh1.*sh2,2));
C2=sum(w.*(sh1(:,2).*sh2(:,1)-sh1(:,1).*sh2(:,2)));

coefmat=[X2,-Y2,W,0;Y2,X2,0,W;Z,0,X2,Y2;0,Z,-Y2,X2];
sol=inv(coefmat)*[X1,Y1,C1,C2]';
a=sol(1:2)';
t=sol(3:4)';
