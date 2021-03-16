function [ dnamask ] = findDnaMask(improt, immask,param)
%Find the largest hole in cell by subtracting intensity values and then
%uses this hole as a seed for adaptive thresholding
%
% Assumes improt and immask (if present) are 3D

%gj july 2012
%gj Oct 9th, 2012  - Uses the distance transform to find nuclear objects
%gj 3/12/13 - comments for clarity. Slight performance enhancements
%gj 8/13/13 - prevent infinite looping in the active contour stage
%rfm 3/9/17 - change call to pc12_cropImg to cropImg
%rfm 3/11/17 - fix masking (did not include third dim in cropping and 
%              restoration) and explicitly test for 3D

% improt = double(improt);
% imsize = size(improt);
% 
% %D. Sullivan check if resolution has been passed in
% %if not, warn user
% if ~isfield(param.model,'resolution')
%     warning('No image resolution given to nuclear hole finding, default scaling 100pixels across');
%     imscale = 100/min(imsize(1:2));
% else
%     %is imscale set?
%     %if not, set it
%     if ~isfield(param,'nucsegresolution')
%         Warning('No nuclear resolution given to nuclear hole finding, default scaling 100pixels across');
%         imscale(1:2) = 100/min(imsize(1:2));
%     else
%         %only resizing in x and y
%         %need to multiply by the downsampling since the images haven't
%         %actually been resized yet, but the resolution has been adjusted
%         imscale = param.model.resolution(1:2)./(param.nucsegresolution(1:2).*param.downsample(1:2));
%     end
% end
% 
% imsizeout = ceil(imsize(1:2).*imscale);
% %D. Sullivan 2/21/13 made this blurring optional
% if ~isfield(param,'holeblur')|| param.holeblur
%     
%     fsize = round(min(imsize(1:2)) / 20);
%     fstd = round(fsize / 5);
%     
%     f = fspecial('gaussian', fsize, fstd);
%     %D. Sullivan 3/20/13 added 'replicate' to avoid edge effects
%     %improt = imfilter(improt,f) +1;
%     improt = imfilter(improt,f,'replicate') +1;
% end
% 
% improtres = imresize(improt, imsizeout);
% 
% if exist('immask', 'var') && ~isempty(immask) && ndims(immask)==2
%     if size(immask,3) == 1;
%         immask = repmat(immask,[1,1,size(improt,3)]);
%     end
%     immaskres = logical(imresize(immask, imsizeout));
%     
%     improtres(~immaskres) = 0;
% elseif exist('immask', 'var') && ~isempty(immask) && all(size(improt)==size(immask))
%     immaskres = logical(imresize(immask, imsizeout));
%     
%     improtres(~immaskres) = 0;
% elseif exist('immask', 'var') && ~isempty(immask) && ~all(size(improt)==size(immask)) 
%     warning(['Mask and protein image are not the same size for the current file. Ignoring mask'])
% % else
% %     %D. Sullivan 3/20/13 need to at least eliminate edge effects of blur.
% %     immaskres = ones(size(improtres)-fstd);
% %     immaskres = padarray(immaskres,[0.5*fstd,0.5*fstd,0.5*fstd]);
% %     improtres(~immaskres) = 0;
% end

if(length(size(improt))<3)
    disp('Warning: image must be 3D in findDnaMask');
    exit
end
improtres = improt;
improtres(~immask) = 0;

[cropimg, range] = cropImg(improtres);
thresh = ml_rcthreshold(cropimg(cropimg > 0));
imthresh = cropimg > thresh;

%find the largest object after thresholding
objs = ml_findobjs(imthresh);
[~, ind] = max(cellfun(@(x) size(x, 1), objs));
obj = double(objs{ind});
ind = sub2ind(size(imthresh), obj(:,1), obj(:,2), obj(:,3));

imthresh = zeros(size(imthresh));
imthresh(ind) = 1;



%D. Sullivan 3/12/13
%Trying to dialate above intensity elements until they enclose the nucleus
% maskimg = cropimg.*imthresh;
% % debug = 1;
% % if debug==1
% %     %show near middle slice b/c it should contain the nucleus
% %     zdisp = floor(size(maskimg,3)/2);
% %     figure, imshow(maskimg(:,:,zdisp),[]);
% % end
% 
% % for i = 1:10
% i = 10;
%     se = strel('ball',i,i);
%     dilatedI = imdilate(maskimg,se);
%     threshD = ml_rcthreshold(dilatedI(dilatedI > 0));
%     imthreshD = dilatedI > threshD;
%     imthresh = imthreshD;
% 
% %     if debug ==1
% %         figure, imshow(imthreshD(:,:,zdisp),[]), figure, imshow(dilatedI(:,:,zdisp),[])
% %         pause
% %     end
%     
% % end

% H = fspecial('disk',10);
% imthresh = imfilter(imthresh,H,'replicate');

for i = 1:size(cropimg,3)
    imthresh(:,:,i) = bwconvhull(imthresh(:,:,i));
end

if exist('param.display', 'var')
    if param.display
        imshow(imthresh(:,:,round(end/2)));
    end
end

if exist('immask', 'var')
    immaskres = immask(range(1):range(2), range(3):range(4), range(5):range(6));
    try 
        imthresh = imthresh .* immaskres;
    catch
        disp('Warning: masking failure in findDnaMask');
        size(imthresh)
        size(immaskres)
    end
end
    
imthresh = cumsum(imthresh,3)>0 & flipdim(cumsum(flipdim(imthresh,3),3),3);

% pad = zeros(size(imthresh));
% pad(2:end-1,2:end-1, 2:end-1) = 1;
% 
% imthresh = imthresh .* pad;

imthresh = padarray(imthresh, [1,1,1]);
cropimg = padarray(cropimg, [1,1,1]);
mask = padarray(immaskres, [1,1,1]);


dist = bwdist(~imthresh | logical(cropimg > thresh) | ~mask);

 t = ml_rcthreshold(dist(logical(imthresh)));

%D. Sullivan 3/20/13 added a param for stiffness 
% dnamask =  region_seg(sqrt(dist), dist > t, 2000,3, 0,0);
dnamask = zeros(size(dist));
if exist('param', 'var') & isfield(param,'stiffness')
    dnamask =  region_seg(sqrt(dist), dist > t, 2000,param.stiffness, ...
        param.display,0.00001,param);
else
    param.stiffness = 3;
    
    maxiter = 5;
    c = 0;
    
    %This will sometimes loop forever if theres no intensity - grj 8/13/13
    while sum(dnamask(:))==0
        c = c+1;
        if c > 5;
            break
        end
        
%         dnamask =  region_seg(sqrt(dist), dist > t, 2000,param.stiffness, 0,0);
        dnamask =  region_seg(sqrt(dist), dist > t, 2000,param.stiffness, ...
            param.display,0.00001,param);
        param.stiffness = param.stiffness/2;

    end
end

% save('dnamasktmp.mat');

%remove padding
dnamask = dnamask(2:end-1, 2:end-1, 2:end-1);

%find the largest object
maskObj = ml_findobjs(dnamask);

[~, ind] = max(cellfun(@(x) size(x,1), maskObj));
maskObj = maskObj{ind};

vox = double(maskObj);
ind = sub2ind(size(dnamask), vox(:,1), vox(:,2), vox(:,3));
dnamask = zeros(size(dnamask));
dnamask(ind) = true;



resmask = zeros(size(improtres));
resmask(range(1):range(2),range(3):range(4),range(5):range(6)) = dnamask;
dnamask = resmask;

%dnamask = imresize(resmask, [size(improt,1), size(improt,2)],'nearest', 'Antialiasing', true);
%dnamask(dnamask > 0) = 1;
%dnamask(dnamask < 0) = 0;

end

