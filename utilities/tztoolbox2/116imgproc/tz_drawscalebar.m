function img2 = tz_drawscalebar(img,res,param)
%TZ_DRAWSCALEBAR Draw a scale bar in an image.
%   IMG2 = TZ_DRAWSCALEBAR(IMG,RES) draws a 10um scale bar at the bottom of
%   the [image] IMG with RES um for each pixel.
%   
%   IMG2 = TZ_DRAWSCALEBAR(IMG,RES,PARAM) customize the scale bar by
%   specifying fields in PARAM:
%       'length' - length of the scale bar. The unit is 'um'. Default: 10um
%       'width'- width of the scale bar. The unit is number of pixels.
%           Default: 10
%       'pos' - position of the scale bar. It is a [point] which specifies
%           the coordinate of the topleft corner of the bar relative to the
%           bottomleft corner of IMG. Default: [20 20].
%   
%   See also

%   08-Feb-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 2
    error('2 or 3 arguments are required')
end

if ~exist('param','var')
    param = struct([]);
end

imageSize = size(img);

param = ml_initparam(param, ...
    struct('length',10,'width',10,'pos',[20,20]));

pixelNumber = round(param.length/res);

maxIntensity = max(img(:));
barpos = [imageSize(1)-param.pos(1), param.pos(2)];
for i=1:size(img,3)
    img(barpos(1):barpos(1)+param.width, ...
        barpos(2):barpos(2)+pixelNumber,i) = maxIntensity;
end

img2 = img;


