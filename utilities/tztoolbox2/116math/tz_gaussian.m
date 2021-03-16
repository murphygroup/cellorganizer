function ydata = tz_gaussian(param,xdata)
%TZ_GAUSSIAN Guassian function for fitting.
%   YDATA = TZ_GAUSSIAN(PARAM,XDATA)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function ydata = tz_gaussian(param,xdata)
%
%OVERVIEW:
%   gaussian functioin for fitting
%PARAMETERS:
%   param - gaussin parameters
%   xdata - data
%RETURN:
%   ydata - y value
%DESCRIPTION
%   f(x)=a*e^(-(x-b)^2/c^2)+d
%HISTORY:
%   ??-???-???? Initial write TINGZ

%val = optimget(my_options,'Display')
ydata = param(1)+param(2)*exp(-(xdata-param(3)).^2./param(4))+param(5)*exp(-(xdata-param(6)).^2./(param(7)));


