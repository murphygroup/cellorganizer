function ydata=tz_projball(param,xdata)
%TZ_PROJBALL Projection of a ball.
%   YDATA = TZ_PROJBALL(PARAM,XDATA)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function ydata=tz_projball(param,xdata)
%
%OVERVIEW:
%   get sum projection of a ball
%PARAMETERS:
%   param - parameters
%   xdata - input coordinates
%RETURN:
%   ydata - output
%DESCRIPTION
%   y=a+b*sqrt(c^2-x^2)
%
%HISTORY:
%   ??-???-2004 Initial write
%   02-NOV-2004 Modified TINGZ

%tz- 27-Sep-2006
% ydata=param(1)+param(2)*sqrt(param(3)-xdata.^2).^(1/2);
%tz-- 

%tz+ 27-Sep-2006
%ydata=param(1)+param(2)*sqrt(1-xdata.^2);
%tz++ 

%tz+ 18-Oct-2006
ydata=param(1)+param(2)*sqrt(1-xdata.^2)+param(3)*sqrt(1-xdata.^2).^(1/2);
%tz++ 
