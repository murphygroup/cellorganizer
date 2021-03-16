function shapemodel = tz_trainshapemodel_ac(shape)
%TZ_TRAINSHAPEMODEL_AC Active shape model.
%   SHAPEMODEL = TZ_TRAINSHAPEMODEL_AC(SHAPE) returns a structure of
%   statistical shape model from the active shape description SHAPE.
%   
%   See also

%   26-May-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University
   
if nargin < 1
    error('Exactly 1 argument is required')
end

shapemodel.sigmas=sqrt(var(shape.pcacom));
shapemodel.covmat=cov(shape.pcacom);
shapemodel.avgshape=shape.avgshape;
shapemodel.pcvec=shape.pcvec;
