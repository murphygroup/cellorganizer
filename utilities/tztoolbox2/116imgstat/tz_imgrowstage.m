function tz_imgrowstage(img,fieldmeth,savedir)
%TZ_IMGROWSTAGE Learning growing stages of an image.
%   TZ_IMGROWSTAGE(IMG,FIELDMETH,SAVEDIR)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function tz_imgrowstage(img,fieldmeth,savedir)
%OVERVIEW
%   learning growing stage of an image
%PARAMETERS
%   img - input image
%   fieldmeth - method of generating a field of pre-stage image
%   savedir - directory for saving results
%RETURN
%   
%DESCRIPTION
%   
%HISTORY
%   28-Apr-2005 Initial write TINGZ
%SEE ALSO
%   tz_trainlogitgrowtex

s=unix(['rm ' savedir '/growstage*.mat']);
s=unix(['mkdir -p ' savedir]);

niter=max(img(:));

for i=1:niter
    growstageimg=img>=niter-i+1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if ~isempty(fieldmeth{1})
%         switch fieldmeth{1}
%         case 'pm4'
%             growstageimg2=bwperim(growstageimg,4);
%         case 'dsm'
%             growstageimg2=bwdist(growstageimg);
%         end
%     else
%         growstageimg2=growstageimg;
%     end
%     
%     if ~isempty(fieldmeth{2})
%         growstageimg2=conv2(growstageimg2,fieldmeth{2},'same');
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    growstageimg=tz_bwproc(growstageimg,fieldmeth{:});
    
%     if any(growstageimg2(:)~=growstageimg(:))
%         keyboard;
%     end
%     
%     imshow(img.*(img>=niter-i+1),[])
    rgbimg=tz_synrgbimg(img>=niter-i,...
        img>=niter-i+1,img>=niter-i+1,'usc');
    imshow(rgbimg,[]);
    
    drawnow
    save([savedir '/growstage' num2str(i) '.mat'],'growstageimg'); 
end