function hfeatures = tz_mtexture(imgs)
%TZ_MTEXTURE Haralick texture of combined images.
%   HFEATURES = TZ_MTEXTURE(IMGS) returns 13 Haralick texture features for
%   the combined image of images in the cell array IMGS.

%   17-May-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argument is required')
end

I=[];
for i=1:length(imgs)
    I=[I,imgs{i}];
end

hfeatures=ml_texture(I); 
hfeatures= hfeatures(1:13,5)';