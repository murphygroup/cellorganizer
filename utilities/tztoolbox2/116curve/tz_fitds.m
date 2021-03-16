function [para,resid]=tz_fitds(ds)
%TZ_FITDS Fit 2D parabola.
%   PARA = TZ_FITDS(DS) returns the parameters of the fitted curve
%   of the targets DS, which is a vector. The value of X coordinate will
%   be from 0 o 1. See TZ_PARABOLA for more details about the curve and
%   parameters.
%   
%   [PARA,RESID] = TZ_FITDS(...) also returns residuals.
%   
%   See also TZ_FITAXLN

%   ??-???-???? Initial write T. Zhao
%   04-NOV-2004 Modified T. Zhao
%       - add comments
%   Copyright (c) Murphy Lab, Carnegie Mellon University

xdata=(0:length(ds)-1)/(length(ds)-1)*2-1;
ydata=ds;

% [maxy,maxy_x]=max(ydata);
% 
% xdata=xdata-max(xdata)/2;%xdata(maxy_x);


[para,resid]=lsqcurvefit(@tz_parabola,[1,63,1,-3,1],xdata,ydata,[1,0,0,-8,0],[3,Inf,Inf,8,4]);

for i=-2:3
    [tmppara,tmpresid]=lsqcurvefit(@tz_parabola,[1,63,1,i,1],xdata,ydata,[1,0,0,-8,0],[3,Inf,Inf,8,4]);
    if tmpresid<resid
        para=tmppara;
        resid=tmpresid;
    end
end
%[para,resid]=lsqcurvefit(@tz_gaussian,[maxy,1,0,1,maxy,1,0],xdata,ydata,[0,0,-Inf,-Inf,0,0,-Inf]);
plot(xdata,ydata,'x');
hold on
plot(xdata,tz_parabola(para,xdata))
hold off

resid=abs(ydata-tz_parabola(para,xdata));