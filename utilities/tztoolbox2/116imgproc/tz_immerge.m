function mergeimg = tz_immerge(img1,img2,method,offset,bg)
%TZ_IMMERGE Merge two images.
%   MERGEIMG = TZ_IMMERGE(IMG1,IMG2,METHOD) returns the merged image of 
%   IMG1 and IMG2. Method specifies how the two images are merged:
%       'max' - maximum value
%       'mean' - mean value
%       'min' - minimal value
%       'mos' - replaced value (replace img2 by img1)
%       'sum' - IMG1+IMG2
%       'sub' - IMG1-IMG2
%   
%   MERGEIMG = TZ_IMMERGE(IMG1,IMG2,METHOD,OFFSET) sets the offset of IMG2.
%
%   MERGEIMG = TZ_IMMERGE(IMG1,IMG2,METHOD,OFFSET,BG) sets the background. The
%   default value is 0;

%   13-Sep-2005 Initial write T. Zhao

if nargin < 3
    error('At least 3 arguments are required')
end

if nargin < 4
    offset = [0;0];
end

if ~exist('bg','var')
    bg = 0;
end

img1 = squeeze(img1);
img2 = squeeze(img2);
mask = ones(size(img1));

img2 = tz_imframe(img2,offset.*(offset>0),[-1,-1],bg);
img1 = tz_imframe(img1,-offset.*(offset<0),[-1,-1],bg);
mask = tz_imframe(mask,-offset.*(offset<0),[-1,-1]);

mergeImgSize = max([size(img1);size(img2)]);

img2 = tz_imframe(img2,[0,0],mergeImgSize(1:2),bg);
img1 = tz_imframe(img1,[0,0],mergeImgSize(1:2),bg);
mask = tz_imframe(mask,[0,0],mergeImgSize(1:2));

switch(method)
    case 'max'
        mergeimg = reshape(max([img1(:),img2(:)],[],2),mergeImgSize);
    case 'min'
        mergeimg = reshape(min([img1(:),img2(:)],[],2),mergeImgSize);
    case 'mean'
        mergeimg = reshape(mean([img1(:),img2(:)],2),mergeImgSize);
    case 'mos'
        mergeimg = img2;
        mergeimg(mask==1) = img1(mask==1);
    case 'sum'
        mergeimg = double(img1)+double(img2);
    case 'sub'
        mergeimg = double(img1)-double(img2);
end
