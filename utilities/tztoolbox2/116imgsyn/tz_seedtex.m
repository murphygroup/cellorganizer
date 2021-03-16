function img = tz_seedtex(imgsize,mask,pkernel,niter)
%TZ_SEEDTEX Use seeding algorithm to generate texture.
%   IMG = TZ_SEEDTEX(IMGSIZE,MASK,PKERNEL,NITER)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function img = tz_seedtex(imgsize,mask,pkernel,niter)
%OVERVIEW
%   Use seeding algorithm to generate texture
%PARAMETERS
%   imgsize - size of the generated image
%   mask -  region where bright pixels could be generated
%   pkernel - transition probabilties
%   niter - number of iteration
%RETURN
%   img - Generated image
%DESCRIPTION
%   
%HISTORY
%   30-Mar-2005 Initial write TINGZ
%SEE ALSO
%   

if isempty(mask)
    mask=ones(imgsize);
end

img=zeros(imgsize);

prevsimg=img;
nbimg=img;

nbkernel=[1 1 1;1 0 1;1 1 1];
pkernel(isnan(pkernel))=0;

for k=1:niter
    k
    updatepos=find(img==0 & mask==1);
    if isempty(updatepos)
        break;
    end
    addpos=find(img>=1 & mask==1);
    nones=nbimg(updatepos);
    curpkernel=pkernel(k,:);
    
    img(updatepos)=binornd(1,curpkernel(nones+1));
    
    if ~isempty(addpos)
        img(addpos)=img(addpos)+1;
    end

    prevsimg=img;
    nbimg=conv2(prevsimg>0,nbkernel,'same');
    
%     img2=zeros(size(img));
%     img2(updatepos)=1;
    imshow(img,[]);
    drawnow
end

% img(1,:)=[];
% img(end,:)=[];
% img(:,1)=[];
% img(:,end)=[];