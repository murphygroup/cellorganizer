function y = ml_frs(x,param)
%ML_FRS Fourier series summation.
%   Y = ML_FRS(X) returns the sum of fouries series to the 5th order,
%   which is f(x)=cos(x)+...+cos(5*x)+sin(x)+...+(5*x).
%   
%   Y = ML_FRS(X,PARAM) specifies how to calculate the sum by specifying
%   the fields in the structure PARAM:
%       'a0' - the offset coefficient. The default value is 0.
%       'order' - the number of components. The default value is 5.
%       'sincoeff' - coeffecients for sin components. It is a vector with
%           the length PARAM.order. The default values are all 1s.
%       'coscoeff' - coeffecients for cos components. It is a vector with
%           the length PARAM.order. The default values are all 1s.
%
%   See also

%   04-Oct-2006 Initial write T. Zhao
%   Copyright (c) 2006 Murphy Lab
%   Carnegie Mellon University
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation; either version 2 of the License,
%   or (at your option) any later version.
%   
%   This program is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   General Public License for more details.
%   
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
%   02110-1301, USA.
%   
%   For additional information visit http://murphylab.web.cmu.edu or
%   send email to murphy@cmu.edu


if nargin < 1
    error('1 or 2 arguments are required');
end

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param,struct('a0',0,'order',10));
param = ml_initparam(param, ...
    struct('sincoeff',ones(1,param.order),'coscoeff',ones(1,param.order)));

y = param.a0;

for n=1:param.order
    y = y+param.sincoeff(n)*sin(n*x);
    y = y+param.coscoeff(n)*cos(n*x);
end


