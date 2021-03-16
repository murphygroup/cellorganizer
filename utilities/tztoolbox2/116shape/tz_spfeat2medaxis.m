function [medaxis,width] = tz_spfeat2medaxis(feat)
%TZ_SPFEAT2MEDAXIS Convert spline features into medial axis.
%   MEDAXIS = TZ_SPFEAT2MEDAXIS(FEAT) returns the medial axis of the shape
%   with medial axis spline features FEAT.
%   
%   [MEDAXIS,WIDTH] = TZ_SPFEAT2MEDAXIS(...) also returns the width.
%   
%   See also

%   28-Dec-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('Exactly 1 argument is required')
end

sp = tz_feats2sp(feat{2},feat{3});

len=round(feat{1});
x=(0:len-1)/(len-1);
axln=spval(sp,x);

sp = tz_feats2sp(feat{4},feat{5});
ds=spval(sp,x);


plot(1:len,axln-ds/2,'.');
hold on
plot(1:len,axln+ds/2,'.');
 
plot(1:len,axln,'x');
hold off
 
axis('equal')
 
% x=1:len;
% figure
% axln2=axln+ds/2;
% gx=[x,x(length(x):-1:1),x(1)];
% gy=[axln-ds/2,axln2(length(x):-1:1),axln(1)-ds(1)/2];
% plot(gx,gy);
% axis([min(x)-10,max(x)+10,min(axln-ds/2)-10,min(axln-ds/2)+max(x)+10])
% axis('equal')
% 
% shape=struct('axln',axln,'ds',ds,'gx',gx,'gy',gy);