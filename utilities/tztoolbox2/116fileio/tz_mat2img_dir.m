function tz_mat2img_dir(res,des,imgtype)
%TZ_MAT2IMG_DIR Convert mat files into image files.
%   TZ_MAT2IMG_DIR(RES,DES) converts mat files under directory RES into
%   BMP files and saves them into directory DES.
%   
%   TZ_MAT2IMG_DIR(RES,DES,IMGTYPE) converts mat files into img files
%   with specified image type, IMGTYPE.
%
%   See also TZ_MAT2IMG_DS

%   20-MAY-2003 Initial write T. Zhao
%   07-JUL-2003 Modified T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if ~exist('imgtype','var')
    imgtype='bmp'
end

filelist=ml_dir([res '/*.mat']);
imgnum=length(filelist);

for i=1:imgnum
    spos=findstr(filelist{i},'/');
    filelist{i}=filelist{i}(max(spos)+1:end);
    filelist{i}
end


for i=1:imgnum
    filename=filelist{i};
    desfilename=[filename(1:end-3) imgtype];
    despath=[des '/' desfilename];
    
    if(~exist(despath))
        img=load([res  '/' filename]);
        grayimg=mat2gray(double(img.proj_image));
        imwrite(grayimg,despath,imgtype);
    end
end
