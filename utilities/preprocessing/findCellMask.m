function [ maskFull] = findCellMask( cellimg, imscale )
%Finds the boundary of a cell in an image of a single cell;
%
%Gregory Johnson, Oct 29, 2012

imsize = size(cellimg);

if ~exist('imscale', 'var')
    imscale = 0.5;
end

%Resize the image to run faster
imres = imresize(cellimg, imscale);

ressize = size(imres);

thresh = ml_rcthreshold(imres(:));

mask = imres > thresh;


[cropimg, range] = pc12_cropImg(imres);

%Get a new mask;
mask = cropimg > thresh;

%find the largest object after thresholding
objs = ml_findobjs(mask);
[~, ind] = max(cellfun(@(x) size(x, 1), objs));
obj = double(objs{ind});
ind = sub2ind(size(mask), obj(:,1), obj(:,2), obj(:,3));

mask = zeros(size(mask));
mask(ind) = 1;


for i = 1:size(mask,3)
    mask(:,:,i) = bwconvhull(mask(:,:,i));
end

imscaled = cropimg;

imscaled(imscaled > thresh) = thresh;

cropmask = region_seg(imscaled, mask, 3000, 0.7, 0, 0);

%Find the largest object
objs = ml_3dfindobj(cropmask);
objs = [objs{:}];
[~,b] = max([objs.size]);

objs = objs(b);

vox = double([objs.voxels]);
voxind = sub2ind(size(cropmask), vox(1,:),vox(2,:), vox(3,:));

cropmask = boolean(zeros(size(cropmask)));
cropmask(voxind) = 1;

%Fill in any holes we find
for i= 1:size(cropmask,3)
    cropmask(:,:,i) = imfill(cropmask(:,:,i), 'holes');
end

cellmask = zeros(size(imres));
cellmask(range(1):range(2), range(3):range(4),:) = cropmask;

maskFull = imresize(cellmask, [size(cellimg,1), size(cellimg,2)], 'Method', 'nearest');
maskFull(maskFull>0) = 1;
maskFull(maskFull<0) = 0;
end

