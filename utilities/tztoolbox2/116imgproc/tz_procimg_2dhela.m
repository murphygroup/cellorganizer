function tz_procimg_2dhela(imgdir,desdir)
%TZ_PROCIMG_2DHELA Preprocess all 2D HeLa images.
%   TZ_PROCIMG_2DHELA(IMAGEDIR,DESDIR)
%   
%   See also

%   17-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function tz_procimg_2dhela(imgdir,desdir)
%OVERVIEW
%   
%PARAMETERS
%   imgdir - 
%   desdir - 
%RETURN
%
%DESCRIPTION
%   
%HISTORY
%   01-Jul-2005 Initial write TINGZ
%SEE ALSO
%   

classes=tz_cleandirs(ml_dir(imgdir));

for i=1:length(classes)
    prot=classes{i};
    protdir=[imgdir '/' prot];
    desprotdir=[desdir '/' prot];
    protfiles= ml_dir([protdir '/prot/*.dat']);
    if ~exist(desprotdir,'dir')
        [status,msg]=mkdir(desdir,prot);
    end
    for j=1:length(protfiles)
        procimgfilename=tz_getfilename_crop_2dhela(protfiles{j},'norm','tif');
        if ~exist([desprotdir '/' procimgfilename],'file')
            protimg=ml_readimage([protdir '/prot/' protfiles{j}]);
            cropfiles{j,1}=tz_getfilename_crop_2dhela(protfiles{j});
            if ~exist([protdir '/crop/' cropfiles{j}],'file')
                cropfiles{j,1}=tz_getfilename_crop_2dhela(protfiles{j},'crop','tiff');
            end
            cropimg=ml_readimage([protdir '/crop/' cropfiles{j}]);
            %     protimg(cropimg==0)=0;
            procimg=ml_preprocess(protimg,cropimg,'ml','yesbgsub');
            normimg=mat2gray(ml_normalize(procimg))*255;
            
            
            imwrite(uint8(normimg),[desprotdir '/' procimgfilename],'tif','Compression','none');
        end
%         imshow(normimg,[]);
%         title(num2str(i));
%         drawnow
    end
end

% title('done')
