function shapemodel = tz_trainshapemodel_sp(spmodels)
%TZ_TRAINSHAPEMODEL_SP Spline shape model (Obsolete).
%   SHAPEMODEL = TZ_TRAINSHAPEMODEL_SP(SPMODELS)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function shapemodel = tz_trainshapemodel_sp(spmodels)
%OVERVIEW
%   
%PARAMETERS
%   spmodels - 
%RETURN
%   shapemodel - 
%DESCRIPTION
%   
%HISTORY
%   26-May-2005 Initial write TINGZ
%SEE ALSO
%   tz_fitmedshape_sp

nshape=length(spmodels);
lens=[];
lnparas=[];
dsparas=[];
for i=1:nshape
    lens=[lens;spmodels{i}.len];
    lnparas=[lnparas;spmodels{i}.lnpara];
    dsparas=[dsparas;spmodels{i}.dspara];
end

shapemodel.axismean=mean([lens,lnparas],1);
shapemodel.axiscov=cov([lens,lnparas]);
shapemodel.dsmean=mean(dsparas,1);
shapemodel.dscov=cov(dsparas);
shapemodel.dssp=spmodels{1}.dssp;
shapemodel.lnsp=spmodels{1}.lnsp;
shapemodel.name='medsp';