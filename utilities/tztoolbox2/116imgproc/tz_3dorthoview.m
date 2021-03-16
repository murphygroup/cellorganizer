function orthimg = tz_3dorthoview(imgs,zmag)
%TZ_3DORTHOVIEW Create orthogonal views for 3D projections.
%   ORTHIMG = TZ_3DORTHOVIEW(IMGS) returns an image which is the 
%   combination of the cell array IMGS. IMGS is a cell array of 3 images,
%   representing Z-projection, X-projection and Y-projection repectively.
%   
%   ORTHIMG = TZ_3DORTHOVIEW(IMGS,ZMAG) also specifies the relative 
%   magnification of Z axis to planar resoution.
%   
%   See also TZ_PROJIMG_3D

%   08-Nov-2005 Initial write T. Zhao
%   Copyright (c) CBI, Carnegie Mellon University

if nargin < 1
    error('1 or 2 arguments are required')
end

if exist('zmag','var')
    imgs{2} = imresize(imgs{2}, ...
        [size(imgs{2},1)*zmag,size(imgs{2},2)],'bilinear');
    imgs{3} = imresize(imgs{3}, ...
        [size(imgs{3},1),size(imgs{3},2)*zmag],'bilinear');
end

gap = round(size(imgs{2},1)/2);

orthimg = zeros(size(imgs{1},1)+gap+size(imgs{2},1), ...
    size(imgs{1},2)+gap+size(imgs{3},2))+double(max(max(imgs{1})));

orthimg(1:size(imgs{1},1),1:size(imgs{1},2)) = imgs{1};
orthimg(size(imgs{1},1)+gap+1:end,1:size(imgs{1},2)) = imgs{2};
orthimg(1:size(imgs{1},1),size(imgs{1},2)+gap+1:end) = imgs{3};