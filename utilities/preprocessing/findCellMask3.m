function [ maskFull] = findCellMask3( cellimg, hessian_sigma )
%Finds the boundary of a cell in an image of a single cell;
%
%Gregory Johnson, Oct 29, 2012
%grj 3/10/14 - returns the largest object

if ~exist('imscale', 'var')
    imscale = 1;
end

if ~exist('hessian_sigma', 'var')
    hessian_sigma = 3;
end

imres = imresize(cellimg, imscale, 'bilinear');

imsize = size(imres);

for i = 1:size(imres,3)
    [dxx, dxy, dyy] = Hessian2D(imres(:,:,i),hessian_sigma);
    [lambda1,lambda2,ix,iy]=eig2image(dxx,dxy,dyy);
    eigim(:,:,i) = max(cat(4, lambda1, lambda2),[],4);
end

imedge = reshape(edge(sqrt(abs(reshape(eigim, [imsize(1), imsize(2)*imsize(3), 1])))), imsize);
    
for i =1:size(imedge,3)
    imedge(:,:,i) = bwmorph(bwmorph(imedge(:,:,i), 'remove'), 'bridge', inf);
    mask(:,:,i) = imerode(imdilate(bwmorph(imedge(:,:,i), 'bridge', inf), strel('disk',7)), strel('disk', 12));
end

mask2 = imfilter(double(mask), fspecial3('gaussian', 7)) > 0.5;


objs = ml_findobjs(mask2);
objsizes = cellfun(@(x) size(x,1), objs);

[~, objind] = max(objsizes);
objs = objs(objind);
% objs = vertcat(objs{:});


cellmask = double(zeros(size(mask2)));
for i = 1:length(objs)
    cropmask = boolean(zeros(size(mask2)));
%     [~,b] = max([objs(i).size]);
% 
%     objs = objs(b);

    vox = double(objs{i}(:,1:3));
    voxind = sub2ind(size(cropmask), vox(:,1),vox(:,2), vox(:,3));


    cropmask(voxind) = 1;


    % cropmask = bwdist(cropmask) <= 3;

    %Fill in any holes we find
    for j= 1:size(cropmask,3)
        cropmask(:,:,j) = imfill(cropmask(:,:,j), 'holes');
    end

    % [~,max_inds,~, min_inds] = extrema(squeeze(sum(sum(cropmask,1),2)));
    % 
    % cell_bottom = 1;
    % if ~isempty(min_inds)
    %     min_inds = sort(min_inds);
    %     cell_bottom = max(min_inds(sort(min_inds) < max_inds(1)));
    % end
    % 
    % cropmask(:,:,1:cell_bottom) = 0;


    cropmask = cumsum(cropmask,3)>0 & flipdim(cumsum(flipdim(cropmask,3),3),3);

    %Fill in any holes we find
    for j = 1:size(cropmask,3)
        cropmask(:,:,j) = imfill(cropmask(:,:,j), 'holes');
    end

    cellmask(cropmask) = i;
end

maskFull = imresize(cellmask, [size(cellimg,1), size(cellimg,2)], 'Method', 'nearest');

corrmat = [];
for i = 1:size(cellimg,3)
    slice = cellimg(:,:,i);
    corrmat(:,i) = slice(:);
end

corrmat = corr(corrmat);
corrinds = find(corrmat(1,:) < 0.995);
maskFull(:,:,1:corrinds(1)-1) = 0;

% maskFull(maskFull>0) = 1;
% maskFull(maskFull<0) = 0;
end

