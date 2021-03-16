function tz_maskdir(dirname,ext,savedir,option)
%TZ_MASKDIR Crop images by mask files and save cropped images.
%   TZ_MASKDIR(DIRNAME,EXT,SAVEDIR,OPTION)
%   
%   See also

%   17-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function mask = tz_maskdir(dirname,option)
%
%OVERVIEW:
%   mask a directory
%PARAMETERS:
%   dirname - image directory for masking
%   ext - extension of the image files
%   savedir - directory for saving masks
%   option - masking for 2d image or 3d image
%RETURN:
%
%DESCRIPTION:
%
%HISTORY:
%   ??-???-???? Initial write TINGZ
%   01-NOV-2004 Modified TINGZ
%       - add comments

switch option
case '3d'
    orgdirname=dirname;
    if dirname(end)=='/'
        dirname(end)=[];
    end
    dirslash=find(dirname=='/');
    if ~isempty(dirslash)
        dirname=dirname(dirslash(end)+1:end);
    end
    if ~exist([savedir,'/' dirname '.mat'],'file')
        img=tz_loadimage(dirname,ext,[]);
        mask=tz_maskimg(img,orgdirname);
        save([savedir,'/' dirname '.mat'],'mask');
    end
case '2d'
    imgfiles=mv_dir([dirname '/*.' ext]);
    for i=1:length(imgfiles)
        img=imread(imgfiles{i});
        mask=tz_maskimg(img);
        save([savedir,'/' imgfiles{i} '.mat'],'mask')
    end
end
    