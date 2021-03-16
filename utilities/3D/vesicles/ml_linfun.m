function y = ml_linfun(x,param)
%ML_LINFUN Linear function. 
%   Y = ML_LINFUN(X,PARAM) returns the linear transform of X. PARAM is a
%   structure with the following parameters:
%       'scale' - scaling
%       'offset' - translation
%   It will give the output Y=X*PARAM.scale+PARAM.offset.  It is recommended
%   to be invertible transformation, which means that PARAM.scale is a
%   scalar or a invertable matrix. If it is not invertable, a psudoinverse
%   matrix will be used.
%   
%   See also ML_GETINVLINFUN

%   10-Aug-2006 Initial write T. Zhao
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


if nargin < 2
    error('Exactly 2 arguments are required');
end

y = x;

if isfield(param,'scale')
    y = y*param.scale;
end

if isfield(param,'offset')
    if length(param.offset)==1
        y = y+param.offset;
    else
        y = ml_addrow(y,param.offset);
    end    
end

