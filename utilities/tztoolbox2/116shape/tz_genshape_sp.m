function shape = tz_genshape_sp(shapemodel)
%TZ_GENSHAPE_SP Generate a shape from medial axis spline model(obsolete).
%   SHAPE = TZ_GENSHAPE_SP(SHAPEMODEL) returns a structure of medial axis
%   shape which has the following fields:
%       'axln' - medial axis
%       'ds' - width
%       'gx' - x coordinates of the boundary
%       'gy' - y coordindates of the boundary
%   
%   See also

%   26-May-2005 Initial write  T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argument is required')
end

lnpara=mvnrnd(shapemodel.axismean,shapemodel.axiscov,1);
len=round(lnpara(1));
x=(0:len-1)/(len-1);
sp=shapemodel.lnsp;
nint=sp.number-sp.order;

if nint>0
    sp.knots(sp.order+1:sp.order+nint)=lnpara(2:nint+1);
end
sp.coefs=lnpara(nint+2:end);
axln=spval(sp,x);

dspara=mvnrnd(shapemodel.dsmean,shapemodel.dscov,1);
nint=sp.number-sp.order;
sp=shapemodel.dssp;
if nint>1
    sp.knots(sp.order+1:sp.order+nint)=dspara(1:nint);
end
sp.coefs=dspara(nint+1:end);
ds=spval(sp,x);


plot(1:len,axln-ds/2,'.');
hold on
plot(1:len,axln+ds/2,'.');
 
plot(1:len,axln,'x');
hold off
 
axis('equal')
 
x=1:len;
figure
axln2=axln+ds/2;
gx=[x,x(length(x):-1:1),x(1)];
gy=[axln-ds/2,axln2(length(x):-1:1),axln(1)-ds(1)/2];
plot(gx,gy);
axis([min(x)-10,max(x)+10,min(axln-ds/2)-10,min(axln-ds/2)+max(x)+10])
axis('equal')

shape=struct('axln',axln,'ds',ds,'gx',gx,'gy',gy);