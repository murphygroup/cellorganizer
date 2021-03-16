function [calpara,xdata,ydata]=tz_texcalrpara(obj,cut)
%TZ_TEXCALRPARA Calibrate the gray level of an object.
%   CALPARA = TZ_TEXCALRPARA(OBJ) returns the calibration parameters of
%   an object to make the object more uniform upon gray levels.
%   
%   CALPARA = TZ_TEXCALRPARA(OBJ,CUT) lets users specify the cut level
%   for removing data, which means pixels with intensities lower than CUT
%   will be considered as background.
%   
%   [CALPARA,XDATA,YDATA] = TZ_TEXCALRPARA(...) also returns data for
%   calibration.

%   06-Sep-2005 Initial write T. Zhao
%   ??-???-???? Initial write T. Zhao
%   05-NOV-2004 Modified T. Zhao
%       - add comments
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('1 or 2 arguments are required')
end

if ~exist('cut','var')
    cut=0.5;
end

img=ml_obj2img(obj,[]);
[calpara,xdata,ydata] = ml_imcalrpara(img,cut);