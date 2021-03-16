function [pkernel,pp,cf] = tz_trainseedtex(img,mask)
%TZ_TRAINSEEDTEX Train seeding texture model.
%   PKERNEL = TZ_TRAINSEEDTEX(IMG,MASK)
%   
%   [PKERNEL,PP,CF] = TZ_TRAINSEEDTEX(...)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function pkernel = tz_trainseedtex(img,mask)
%OVERVIEW
%   
%PARAMETERS
%   img - 
%   mask - 
%RETURN
%   pkernel - 
%DESCRIPTION
%   
%HISTORY
%   30-Mar-2005 Initial write TINGZ
%SEE ALSO
%

if isempty(mask)
    mask=ones(size(img));
end

img=double(img);
imgsize=size(img);
img(find(mask==0))=0;
niter=max(img(:));

pkernel=zeros(1,9);
tkernel=zeros(1,9);
prevsimg=(img==niter);

nbimg=zeros(size(img));
nbkernel=[1 1 1;1 0 1;1 1 1];
nbimg=conv2(prevsimg,nbkernel,'same');
maskpos=find(mask==1);
nmaskpixel=length(maskpos);
for k=niter-1:-1:1
    k

    prevsimg=(img>k);
    cursimg=(img==k);
    prevpkernel=pkernel;
    prevtkernel=tkernel;
    
    curmaskpos=maskpos;
    curmaskpos(prevsimg(maskpos)==1)=[];
    
    nones=nbimg(curmaskpos);
    npost=tz_label2post(nones+1,9);
    tkernel=tkernel+sum(npost,1);
    pkernel=pkernel+sum(npost(cursimg(curmaskpos)==1,:),1);
    
    
    nbimg=conv2(cursimg,nbkernel,'same')+nbimg;
    curpp=(pkernel-prevpkernel)./(tkernel-prevtkernel);
    curpp(tkernel-prevtkernel==0)=NaN;
    pp(k,:)=curpp;
    imshow(prevsimg,[])
    cf(k)=sum(prevsimg(:));
    %     plot(pp(:,3));
    %     drawnow
    [pkernel;tkernel];
end

pkernel=pkernel./tkernel;