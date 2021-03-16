function ml_immesh(img,param)
%ML_IMMESH Show image as a 3-D mesh surface.
%   ML_IMMESH(IMG,PARAM) plots the image IMG as a mesh where (X,Y) are the
%   pixel coordinates and Z is the intensity of the image. PARAM is a structure
%   with the following fields:
%     'step' - step of making a mesh
%     'plot' - plot method. It supports 'mesh','surf','surface', which 
%          correspond to the Matlab functions with the same names.
%     'plotparam' - plot parameters. See the help file of the functions in 
%          Matlab. It is a cell array containing a list of parameters.
%   
%   See also

%   03-Nov-2006 Initial write T. Zhao
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

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param,struct('step',1,'plot','mesh','plotparam',{{}}));

img = img';
[x,y] = meshgrid(1:param.step:size(img,1),1:param.step:size(img,2));

z = img(sub2ind(size(img),x,y));

switch param.plot
    case 'mesh'
        mesh(x,y,z,param.plotparam{:});
    case 'surf'
        surf(x,y,z,param.plotparam{:});
    case 'surface'
        surface(x,y,z,param.plotparam{:});
    otherwise
        error(['Unrecognized plot method: ' param.plot]);
end

