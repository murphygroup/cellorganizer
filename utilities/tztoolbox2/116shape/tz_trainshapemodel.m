function shapemodel = tz_trainshapemodel(shapes,param)
%TZ_TRAINSHAPEMODEL Train statistical shape model from shapes.
%   SHAPEMODEL = TZ_TRAINSHAPEMODEL(SHAPES) returns a statistical meidal 
%   axis shape model from the training data SHAPES, which is a cell array 
%   of shapes of the same format of 'mxp'.
%
%   SHAPEMODEL = TZ_TRAINSHAPEMODEL(SHAPES,PARAM) lets the user specify the
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
%       'isshow' - show intermediate results if it is not 0. Currently is
%           is only useful for 'act'.
%
%   See also TZ_GENSHAPE

%   02-Jan-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('1 or 2 arguments are required')
end
warning(tz_genmsg('of','tz_trainshapemodel','ml_trainshapemodel'));
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
    case 'mxp'
        shapemodel.name = 'mxp';
        [knots1,coef1] = tz_sp2feat(shapes{1}.spmedaxis);
        [knots2,coef2] = tz_sp2feat(shapes{1}.spwidth);
        shapemodel.nmedknots = length(knots1);
        shapemodel.nwidthknots = length(knots2);  
        
%         [knots1,coef1] = tz_sp2feat(shapes{1}.spmedaxis);     
%         [knots2,coef2] = tz_sp2feat(shapes{1}.spwidth);
%         
       
%         medfeats = [shapes{1}.length,knots1,coef1];
%         widthfeats = [knots2,coef2];
        
        for k=1:length(shapes)
            [knots1,coef1] = tz_sp2feat(shapes{k}.spmedaxis);
            [knots2,coef2] = tz_sp2feat(shapes{k}.spwidth);
            if ~isempty(param.constknots{1})
                knots1 = [];
            end
            if ~isempty(param.constknots{2})
                knots2 = [];
            end
            medfeats(k,:) = [shapes{k}.length,knots1,coef1];
            widthfeats(k,:) = [knots2,coef2];
        end
        
        shapemodel.medmean = mean(medfeats,1);
        shapemodel.medcov = cov(medfeats);
        shapemodel.widthmean = mean(widthfeats,1);
        shapemodel.widthcov = cov(widthfeats);
        shapemodel.constknots = param.constknots;
    case 'act'  %active shape model
        param = tz_initparam(param,struct('niter',50,'isshow',0));
        shape = tz_averageshape(shapes,param.niter,param.isshow);
        shapemodel = tz_trainshapemodel_ac(shape)
        shapemodel.name = 'act';
    case 'rds'
        if ~isfield(param,'tz_rdistpca')
            param.tz_rdistpca = struct([]);
        end
        shapemodel.name = 'rds';
        shapemodel.stat = tz_rdistpca(shapes,param.tz_rdistpca);
    otherwise
        shapemodel.name = 'unknown';
end
