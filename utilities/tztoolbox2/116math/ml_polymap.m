function y = ml_polymap(x,param)
%ML_POLYMAP Polynomial calculation
%   Y = ML_POLYMAP(X) returns [1 X].
%   
%   Y = ML_POLYMAP(X,PARAM) returns the polynomial results of X. PARAM 
%   is a structure that determines the polynomial. It has the following
%   fields:
%       'order' - order of the polynomial
%       'const' - constant. If it is empty, no constant will be added. The
%            default value is 1.
%       'coeff' - coeffients (do not include constant). If it is empty, all
%            coffecients are 1s. Otherwise, it must be a row
%            vector and have the length ORDER. The default value is empty.
%   
%   See also

%   04-Sep-2006 Initial write T. Zhao
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
    param = struct('order',1);
end

if ~isfield(param,'order')
    error('The order is not specified');
end

param = ml_initparam(param,struct('const',1, ...
    'coeff',ones(1,param.order)));

y = [];

for i=1:param.order
    xx = x.^i;
    if ~isempty(param.coeff)
        xx = x.^i*param.coeff(i);
    end
    
    y = [y xx];
end

if ~isempty(param.const)
    y = [zeros(size(x,1),1)+param.const,y];
end






