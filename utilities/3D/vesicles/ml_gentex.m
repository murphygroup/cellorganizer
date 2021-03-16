function img = ml_gentex(model,param)
%ML_GENTEX Texture synthesis.
%   IMG = ML_GENTEX(MODEL) returns a texture image that is synthesized
%   from the texture model MODEL, which a structure returned from
%   ML_TEXMODEL. The image size is [256 256].
%   
%   IMG = ML_GENTEX(MODEL,PARAM) allows customizing the texture synthesis
%   by specifying the structure PARAM with the following fields:
%       'size' - texture size [row column]
%   Some other fields depends on MODEL.name:
%       'pyr' : pyramid model
%           'Niter' - number of iteration.
%   
%   See also ML_TEXMODEL

%   20-Aug-2006 Initial write T. Zhao
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

param = ml_initparam(param,struct('size',[256 256]));

switch model.name
    case 'pyr'
        param = ml_initparam(param,struct('Niter',25));
        img = textureSynthesis(model,param.size,param.Niter);
    otherwise
        error('Unrcognized texture model');
end
