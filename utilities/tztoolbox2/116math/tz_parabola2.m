function ydata = tz_parabola2(param,xdata)
%TZ_PARABOLA2 Ploynomial function for fitting.
%   YDATA = TZ_PARABOLA2(PARAM,XDATA)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University


%function tz_parabola2(param,data)
%
%OVERVIEW:
%   polynomial
%ydata=sqrt(param(2)-param(1)*xdata.^2)+param(3)
ydata = param(5)*xdata.^4+param(1)*xdata.^3+param(2)*xdata.^2+param(3).^xdata+param(4);
%ydata=param(1)+param(2)*gamma(param(3)*xdata+param(4));param(1)+
%param(1)*xdata.^4+ +++param(2)*xdata+param(3)*xdata.^2