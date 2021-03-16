function img2 = tz_normreg(img)
%TZ_NORMREG Convert a [region image] to a [normalized region image].
%   IMG2 = TZ_NORMREG(IMG) returns a [normalized region image] which
%   has the same regions as those in the [region image] IMG.
%   
%   See also

%   23-Jan-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('Exactly 1 argument is required')
end

bwimg = img>0;
newLabeledImage = bwlabel(bwimg);

newLabeledImage = immultiply(newLabeledImage, ...
    zeros( size(newLabeledImage) )+double( max( img(:) ) ) );

img2 = imadd(double(img),newLabeledImage);

eval(['img2 = ' class(img) '(img2);']);

img2(img2>0) = tz_squeezenum(img2(img2>0));
