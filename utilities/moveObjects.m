function newpos = moveObjects(pos,Dc,dt,param)
%MOVEOBJECTS This function moves object positions according to the specified parameters
%
%Inputs:
%pos = nx3 array of current positions for n objects in 3D
%Dc = diffusion coefficient (um^2/second)
%dt = time step
%
%Outputs:
%
%Written by: Devin Sullivan 7/18/13

% Copyright (C) 2013 Murphy Lab
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

param = ml_initparam(param,struct('method','brownian'));

switch param.method
    case 'brownian'
        newpos = pos+randn(size(pos,1),3)*Dc*sqrt(dt);
        
    otherwise
        warning('Unrecognized object movement method, returning initial positions.');
        newpos = pos;
end