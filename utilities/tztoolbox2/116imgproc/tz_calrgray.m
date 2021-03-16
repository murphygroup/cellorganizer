function img2=tz_calrgray(img,para,isinv)
%TZ_CALRGRAY Calibrate gray level of pixels according to positions.
%   IMG2 = TZ_CALRGRAY(IMG,PARA,ISINV) returns an image that is the
%   calibration of the [image] IMG. The calibartion is to adjust the image
%   intensity according to the postion relative to edge. The paramters of
%   calibration is PARA, which is a row vection with length 3. See
%   TZ_PROJBALL for more details.If ISINV is set to 1, the function will
%   try the inverse procedure. In IMG the pixels with intensity no greater
%   than 0 will be taken as background.
%   
%   See also

%   ??-???-2004 Initial write TINGZ
%   03-NOV-2004 Modified TINGZ
%       - add comments
%   Copyright (c) Murphy Lab, Carnegie Mellon University


img2=tz_edgedist(img);
distvec=tz_mat2vec(img2);
imgvec=tz_mat2vec(img);

distvec=(max(distvec)-distvec)/max(distvec);
egray=tz_projball(para,distvec);

eratio=max(egray)./egray(egray>0);
if isinv==1
    eratio=1./eratio;
end

imgvec(egray>0)=imgvec(egray>0).*eratio;

img2=reshape(imgvec,size(img));
