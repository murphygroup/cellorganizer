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
% Copyright (C) 2007  Murphy Lab
% Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.
%
% For additional information visit http://murphylab.web.cmu.edu or
% send email to murphy@cmu.edu

ydata=param(1)+param(2)*sqrt(1-xdata.^2)+param(3)*sqrt(1-xdata.^2).^(1/2);
%tz++ 
