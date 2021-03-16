function [lnpara,resid,y]=tz_fitaxln(axln)
%TZ_FITAXLN Fit a curve by guassion mixture function.
%   LNPARA = TZ_FITAXLN(AXLN) returns the parameters of the fitted curve
%   of the point array AXLN, which is a Nx2 matrix if there are N points.
%   The 1st column and 2nd column of AXLN are X and Y coordinate
%   repectively. See TZ_GAUSSIAN for more details about the curve and
%   parameters.
%   
%   [LNPARA,RESID,Y] = TZ_FITAXLN(...) also returns the residuals and 
%   fitted values.
%   
%   See also TZ_FITDS

%   ??-???-2004 Initial write T. Zhao
%   04-NOV-2004 Modified T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

xdata=axln(:,1);
ydata=axln(:,2);
ydata=ydata-min(ydata);
xdata=(0:length(xdata)-1)/(length(xdata)-1)*2-1;
xdata=xdata';


lb=[-Inf,0,-Inf,0,0,-Inf,0];
ub=[Inf,Inf,Inf,2,Inf,Inf,2];
[lnp1,lnr1]=lsqcurvefit(@tz_gaussian,[1,1,0,1,0,0,1],xdata,ydata,lb,ub);
ydata2=-ydata;
ydata2=ydata2-min(ydata2);
[lnp2,lnr2]=lsqcurvefit(@tz_gaussian,[1,1,0,1,0,0,1],xdata,ydata2,lb,ub);

if lnr1<lnr2
    lnpara=lnp1;
    resid=lnr1;
    y=ydata;
else
    lnpara=lnp2;
    resid=lnr2;
    y=ydata2;
end

if lnpara(3)>lnpara(6)
    lnpara=[lnpara(1),lnpara(5:7),lnpara(2:4)];
end

if lnpara(2)<1e-5
    lnpara(2:4)=0;
end

if lnpara(5)<1e-5
    lnpara(5:7)=0;
end

plot(axln(:,1),y,'r.');
hold on
plot(axln(:,1),tz_gaussian(lnpara,xdata))
axis([min(axln(:,1)),max(axln(:,1)),-64,64])
hold off

resid=abs(y-tz_gaussian(lnpara,xdata));
resid=resid';