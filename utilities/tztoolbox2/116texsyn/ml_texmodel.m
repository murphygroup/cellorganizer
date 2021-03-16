function tex = ml_texmodel(img,param)
%ML_TEXMODEL Train a texture model.
%   TEX = ML_TEXMODEL(IMG) returns a texture model that is trained on the
%   image IMG. See TEXTUREANALYSIS about more details of the model.
%   
%   TEX = ML_TEXMODEL(IMG,PARAM) allows users to customize the way of
%   training the texture model by PARAM, which is a structure with the
%   following fields: 
%       'modelname' - name of the model
%       'preprocess' - preprocessing function. see ML_IMPROC for details.
%   
%   See also

%   18-Aug-2006 Initial write T. Zhao
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
    param = struct('modelname','pyr','preprocess',[]);
end

param = ml_initparam(param,struct('modelname','pyr'));
if ~isempty(param.preprocess)
    img2 = ml_improc(img,param.preprocess);
end

switch param.modelname
    case 'pyr'
        param = ml_initparam(param,struct('Nsc',4,'Nor',4,'Na',7));
        tex = textureAnalysis(img2,param.Nsc,param.Nor,param.Na);
        tex.name = param.modelname;
    otherwise
        error('Unrecognized texture model');
end
