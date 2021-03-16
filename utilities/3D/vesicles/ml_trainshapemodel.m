function shapemodel = ml_trainshapemodel(shapes,param)
%TZ_TRAINSHAPEMODEL Train statistical shape model from shapes.
%   SHAPEMODEL = ML_TRAINSHAPEMODEL(SHAPES) returns a statistical meidal 
%   axis shape model from the training data SHAPES, which is a cell array 
%   of shapes of the same format of 'mxp'.
%
%   SHAPEMODEL = ML_TRAINSHAPEMODEL(SHAPES,PARAM) lets the user specify the
%   model name and training parameters by PARAMA, which is a structure that
%   could have the following fields:
%       'modelname' - model name
%           'mxp' : medial axis model. SHAPES is a cell array of objects 
%               with medial axis spline representation
%               'constknots' - constant knots
%           'act' : active shape model. SHAPES is a cell array of [point
%               array]s.
%           'rds' : distance ratio shape model. SHAPES is a cell array of
%               structures of cell codes.
%       'niter' - number of iteration (See TZ_AVERAGESHAPE). It is only
%           useful for 'act'. The default value is 50.
%       'isshow' - show intermediate results if it is not 0. Currently it
%           is only useful for 'act'.
%
%   See also ML_GENSHAPE

%   02-Jan-2006 Initial write T. Zhao

% Copyright (C) 2006  Murphy Lab
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

if nargin < 1
    error('1 or 2 arguments are required')
end

if isempty(shapes)
    shapemodel.name = 'unknown';
    return;
end

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param,struct('modelname','mxp', ...
    'constknots',{{[],[]}}));

switch param.modelname
    case 'mxp' %medial axis spline
        shapemodel.name = 'mxp';
        
        %knots and cofficients of the spline
        [knots1,coef1] = ml_sp2feat(shapes{1}.spmedaxis);
        [knots2,coef2] = ml_sp2feat(shapes{1}.spwidth);
        
        shapemodel.medaxis.nknots = length(knots1);
        shapemodel.width.nknots = length(knots2);
        
        shapemodel.medaxis.constknot = param.constknots{1};
        shapemodel.width.constknot = param.constknots{2};
        
%         if ~isempty(param.constknots{1})
%             shapemodel.medaxis.constknot = [];
%         end
%         if ~isempty(param.constknots{2})
%             shapemodel.width.constknot = [];
%         end
        %shapemodel.nmedknots = length(knots1);
        %shapemodel.nwidthknots = length(knots2);  
                
        shapemodel.medaxis.name = 'bspline';
        shapemodel.width.name = 'bspline';
        for k=1:length(shapes)
            [knots1,coef1] = ml_sp2feat(shapes{k}.spmedaxis);
            [knots2,coef2] = ml_sp2feat(shapes{k}.spwidth);

            if isempty(param.constknots{1})
                medfeats(k,:) = [shapes{k}.length,knots1,coef1];
            else
                medfeats(k,:) = [shapes{k}.length,coef1];
            end
            if isempty(param.constknots{2})
                widthfeats(k,:) = [knots2,coef2];
            else
                widthfeats(k,:) = coef2;
            end
        end
        
        shapemodel.medaxis.stat = ml_estpdf(medfeats,struct('name','mvn'));
        shapemodel.width.stat = ml_estpdf(widthfeats,struct('name','mvn'));
%         shapemodel.medmean = mean(medfeats,1);
%         shapemodel.medcov = cov(medfeats);
%         shapemodel.widthmean = mean(widthfeats,1);
%         shapemodel.widthcov = cov(widthfeats);
%         shapemodel.constknots = param.constknots;
    case 'act'  %active shape model
        param = ml_initparam(param,struct('niter',50,'isshow',0));
        shape = ml_averageshape(shapes,param.niter,param.isshow);
        shapemodel = ml_trainshapemodel_ac(shape)
        shapemodel.name = 'act';
    case 'rds'
        if ~isfield(param,'ml_rdistpca')
            param.ml_rdistpca = struct([]);
        end
        shapemodel.name = 'rds';
        shapemodel = ml_initparam(shapemodel, ...
            ml_rdistpca(shapes,param.ml_rdistpca));
    otherwise
        shapemodel.name = 'unknown';
end
