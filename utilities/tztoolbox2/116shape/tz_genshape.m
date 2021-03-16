function shape = tz_genshape(shapemodel,param)
%TZ_GENSHAPE Randomly generates a shape.
%   SHAPE = TZ_GENSHAPE(SHAPEMODEL) returns a shape that is randomly
%   sampled from the statistical shape model SHAPEMODEL.
%   
%   SHAPE = TZ_GENSHAPE(SHAPEMODEL,PARAM) spedifies the parameters for
%   generation. Currently it is only useful for 'act' model. THe field is
%   'ncomp', which is the number of principal components for generation.
%   
%   See also tz_trainshapemodel

%   02-Jan-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('1 or 2 arguments are required')
end

warning(tz_genmsg('of','tz_genshape','ml_genshape'));

if ~exist('param','var')
    param = struct([]);
end

switch shapemodel.name
    case 'mxp'
        shape.format = 'mxs';
        medfeat = mvnrnd(shapemodel.medmean,shapemodel.medcov,1);
        widthfeat = mvnrnd(shapemodel.widthmean,shapemodel.widthcov,1);
        shape.length = round(medfeat(1));
        if isempty(shapemodel.constknots{1})
            knots1 = medfeat(2:shapemodel.nmedknots+1);
            coef1 = medfeat(shapemodel.nmedknots+2:end);
        else
            knots1 = shapemodel.constknots{1};
            coef1 = medfeat(2:end);
        end
        
        if isempty(shapemodel.constknots{2})
            knots2 = widthfeat(1:shapemodel.nwidthknots);
            coef2 = widthfeat(shapemodel.nwidthknots+1:end);
        else
            knots2 = shapemodel.constknots{2};
            coef2 = widthfeat(1:end);
        end
        
        shape.spmedaxis = tz_feat2sp(knots1,coef1);
        shape.spwidth = tz_feat2sp(knots2,coef2);
    case 'act'
        param = tz_initparam(param,struct('ncomp',5));
        shape.format = 'crd';
        ncom = param.ncomp;
        ss = shapemodel.avgshape(:)+shapemodel.pcvec(:,1:ncom)* ...
            mvnrnd(zeros(1,ncom),shapemodel.covmat(1:ncom,1:ncom),1)';
        shape.pts=[ss(1:360),ss(361:end)];
    otherwise
        error('Unrecognized shape model.');
end
