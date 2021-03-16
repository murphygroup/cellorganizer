function img = ...
    tz_growtex(imgsize,beta,option,cellbody,nucbody,niter,seedimg,fieldmeth)
%TZ_GRWOTEX Synthesize image by growing texture model.
%   IMG = 
%   TZ_GRWOTEX(IMGSIZE,BETA,OPTION,CELLBODY,NUCBODY,NITER,SEEDIMG,FIELDMETH)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function img = tz_growtex(imgsize,beta,option,cellbody,nucbody,niter,seedimg,nbkernel)
%OVERVIEW
%   grow texture
%PARAMETERS
%   imgsize - image size
%   beta - logistic model parameters
%   option - option for logistic variable transformation
%   cellbody - cell body
%   nucbody - nucleus
%   niter - number of iterations (growing stages)
%   seedimg - seed
%   nbkernel - kernel for neighbor's influence
%RETURN
%   img - generated image
%DESCRIPTION
%   
%HISTORY
%   13-Apr-2005 Initial write TINGZ
%SEE ALSO
%   
 
if ~isempty(cellbody)    
    cellmask=tz_obj2img(cellbody,imgsize,{'2d','bn'});
    celledgeimg=bwperim(cellmask,8);
    celldistimg=bwdist(celledgeimg);
else
    cellmask=ones(imgsize);
    celldistimg=[];
end

if ~isempty(nucbody)
    dnamask=tz_obj2img(nucbody,imgsize,{'2d','bn'});
    dnaedgeimg=bwperim(dnamask,8);
    dnadistimg=bwdist(dnaedgeimg);
else
    dnamask=zeros(imgsize); 
    dnadistimg=[];
end

mask=cellmask-dnamask;
maxcelldist=max(celldistimg(mask(:)==1));

if isempty(seedimg)
    img=zeros(imgsize);
else
    img=seedimg;
end

prevsimg=img;
nbimg=img;

% if ~exist('nbkernel','var')
%     nbkernel=[];  
% end
% if isempty(nbkernel)
%     nbkernel=[1 1 1;1 0 1;1 1 1];
% end

% for k=1:niter
k=0;
while k<=niter
    k
    prevsimg=img;
    
%     if isempty(nbkernel)
%         nbimg=bwdist(prevsimg>0);
%     else
%         nbimg=conv2(prevsimg>0,nbkernel,'same');
%     end
    
    nbimg=tz_bwproc(prevsimg>0,fieldmeth{:});
    
    updatepos=find(img==0 & mask==1);
    if isempty(updatepos)
        break;
    end
    addpos=find(img>=1 & mask==1);
    nones=nbimg(updatepos);

    if ~isempty(cellbody)
        celldists=celldistimg(updatepos);
    else
        celldists=[];
    end
    if ~isempty(nucbody)
        dnadists=dnadistimg(updatepos);
    else
        dnadists=[];
    end
    
    zs=zeros(length(nones),1);
    %     img(updatepos)=binornd(1,...
    %             tz_evallogistic([nones,celldists,dnadists,niter-k+1+zs,1+zs],beta));
    if strcmp(option,'gt11b')
        celldists=celldists/maxcelldist;
        dnadists=dnadists/maxcelldist;
    end
    x=tz_mapdata([1+zs,k+zs,celldists,dnadists,nones],option);
    img(updatepos)=binornd(1,tz_evallogistic(x,beta));
%    keyboard
    if ~isempty(addpos)
        img(addpos)=img(addpos)+1;
    end
        
    imshow(img,[]);
    drawnow
    k=k+1;
    
    if all(img==0)
        k=0;
    end
end 
