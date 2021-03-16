function teximg = tz_imgproc4tex(img,cropimage,scale,har_intbins,har_pixsize)
%TZ_IMGPROC4TEX Get the image ready for texture feature calculation.
%   TEXIMG = TZ_IMGPROC4TEX(IMG,CROPIMAGE,SCALE,HAR_INTBINS,HAR_PIXSIZE)
%   
%   See also

%   17-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function teximg = tz_imgproc4tex(img,cropimage,scale,har_intbins,har_pixsize)
%OVERVIEW
%   
%PARAMETERS
%   img - original image
%   cropimage - crop mask
%   scale - micrometers per pixel
%   har_intbins - the number of intensity bins to use for Haralick texture features
%   har_pixsize - the pixel size (in micrometers) to use for Haralick texture features
%RETURN
%   teximg - image for texture feature calculation
%DESCRIPTION
%   
%HISTORY
%   17-May-2005 Initial write TINGZ
%SEE ALSO
%   tz_mtexture

bgsub = 'yesbgsub';

if isempty(scale)
    scale=0;
end

way = 'ml';

% Preprocess the image
[procimage, prot_maskimage,nonobjimg] = ml_preprocess( img, cropimage, ...
						  way, bgsub);% Calculate the scale factor
                      
default_scale = 0.23; %%% 0.23um/pixel
if( scale == 0)
    scale = default_scale;
    scale_factor = 1;
else
    scale_factor = scale / default_scale;
end;

DEFAULT_INTENSITY_BINS = 256;
DEFAULT_HAR_PIXEL_SIZE = 1.15;

% Check arguments & use defaults where necessary
if isempty(har_intbins)
    har_intbins = DEFAULT_INTENSITY_BINS;
end
if isempty(har_pixsize)
    har_pixsize = DEFAULT_HAR_PIXEL_SIZE;
end

har_rescale_factor = scale / har_pixsize;
resized_img = imresize(procimage, har_rescale_factor, 'bilinear');
teximg = uint8(ml_imgscale(resized_img)*har_intbins);

% imshow(teximg,[])