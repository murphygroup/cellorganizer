function fi = ml_getinvlinfun(f)
%ML_GETINVLINFUN Get inverse function of a linear function
%   FI = ML_GETINVLINFUN(F) returns the inverse function of F, which must be
%   a ML_LINFUN fucntion.
%    
%   See also ML_LINFUN

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


if nargin < 1
    error('Exactly 1 argument is required');
end

if ~strcmp(f.funname,'ml_linfun')
    error('The input funciton must be a ml_linfun');
end

fi.funname = f.funname;

if isfield(f.param,'scale')
    if size(f.param.scale,1)==1
        fi.param.scale = 1./f.param.scale;
    else
        fi.param.scale = pinv(f.param.scale);
    end
end

if isfield(f.param,'offset')
    if isfield(f.param,'scale')
        if size(f.param.scale,1)==1
            fi.param.offset = -f.param.offset./f.param.scale;
        else
            fi.param.offset = -f.param.offset*pinv(f.param.scale);
        end
    else
        fi.param.offset = -f.param.offset;
    end
end

    
