function tz_imshowdir(imdir,ext,mode,mask)
%TZ_IMSHOWDIR Shows images under a directory.
%   TZ_IMSHOWDIR(IMDIR,EXT,MODE) shows images with the extention EXT under
%   the directory IMDIR. MODE specifies how to show the images:
%       'ID' - show images one by one indepently, i.e. each image is scaled
%              based on itself.
%       '3D' - show images one by one but they are scaled based on all of
%              images.
%       'proj' - show projections of the images.
%   
%   TZ_IMSHOWDIR(IMDIR,EXT,MODE,MASK) specifies a 2D mask for showing.

%   04-MAR-2003 Initial write T. Zhao
%   25-MAR-2003 Modified T. Zhao
%   20-MAY-2003 Modified T. Zhao
%   14-DEC-2003 Modified T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 3
    error('3 or 4 arguments are required')
end

if ~exist('mask','var')
    mask=[];
end

switch mode
    case 'ID',
        [img,imgfiles]=tz_loadimage(imdir,ext,65535);
        no_of_slices=size(img,3);
        
        for slice_no = 1:no_of_slices
            slice_no
            slice = img(:,:,slice_no);
            if ~isempty(mask)
                slice=slice.*mask;
            end
            
            imshow(slice,[min(slice(:)),max(slice(:))]);
            title(imgfiles{slice_no})
            pause
        end
                
    case '3D',
        img=tz_loadimage(imdir,ext,65535);
        minpixel=min(img(:))
        maxpixel=max(img(:))
        no_of_slices=size(img,3);
        [x,y,z]=find(img==maxpixel);
        for slice_no = 1:no_of_slices
            slice_no
            slice = img(:,:,slice_no);
            [x,y]=find(slice==maxpixel);
            if ~isempty(mask)
                slice=slice.*mask;
            end
            imshow(slice,[0,maxpixel]);
            pause
        end
    case 'proj',
        img1=tz_loadimage(imdir,ext,65535);
        img2=double(img1);
        img2proj=sum(img2,3);
        if ~isempty(mask)
            img2proj=img2proj.*mask;
        end
        image(img2proj);
        title(imdir);
        test=1;
    otherwise,
        warning('Invalid mode');
end
