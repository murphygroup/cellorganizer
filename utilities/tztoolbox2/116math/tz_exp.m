function ydata = tz_exp(param,xdata)
%TZ_EXP PDF of exponential distribution.
%   YDATA = TZ_EXP(PARAM,XDATA)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function ydata = tz_exp(param,xdata)
%OVERVIEW
%   exponential function y=a*exp(bx)+c
%PARAMETERS
%   param - parameters
%   xdata - variables
%RETURN
%   ydata - return value
%DESCRIPTION
%   ydata=param(1)*exp(xdata*param(2))+param(3)
%HISTORY
%   04-Apr-2005 Initial write TINGZ
%SEE ALSO
%   

ydata=param(1)*exp(xdata*param(2))+param(3);