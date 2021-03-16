function tz_tif2bmp_dir(res,des,ext,is_3d,imgtype)
%TZ_TIF2BMP_DIR Convert images into another image format.
%   TZ_TIF2BMP_DIR(RES,DES,EXT,IS_3D) converts images under directory RES
%   into BMP images and saves them under directory DES. EXT is the
%   extension of the original images. The images will be scaled
%   automacally. If IS_3D is not 0, the images will be scaled upon all
%   images under the directory.
%   
%   TZ_TIF2BMP_DIR(RES,DES,EXT,IS_3D,IMGTYPE) specifies the converted
%   image format, which is BMP for default.
%   
%   See also TZ_TIF2BMP_DS

%   ??-???-???? Initial write T. Zhao
%   12-JUL-2003 T. Zhao
%   14-DEC-2003 T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 4
    error('4 or 5 arguments are required')
end

if ~exist('imgtype','var')
    imgtype='bmp'
end

filelist=ml_dir([res '/*.' ext]);
imgnum=length(filelist);

if imgnum==0
    return ;
end

for i=1:imgnum
    spos=findstr(filelist{i},'/');
    if ~isempty(spos)
        filelist{i}=filelist{i}(max(spos)+1:end);
    end
    filelist{i}
end

if is_3d==0
    for i=1:imgnum
        filename=filelist{i};
        tifimg=ml_readimage([res  '/' filename]);
        grayimg=mat2gray(double(tifimg));
        desfilename=[filename(1:end-3) imgtype];
        imwrite(grayimg,[des '/' desfilename],imgtype);
    end
else 
        
    img=tz_loadimage(res,ext,65535);
    
    max3d=max(img(:))
    
    %grayimg=mat2gray(double(img));
    for i=1:imgnum
        maxslice=max(max(img(:,:,i)));
        filename=filelist{i};
        slice=img(:,:,i);
        grayimg=mat2gray(double(slice),[0 double(max3d)]);
        desfilename=[filename(1:end-3) imgtype];
        imwrite(grayimg,[des '/' desfilename],imgtype);
    end
end
