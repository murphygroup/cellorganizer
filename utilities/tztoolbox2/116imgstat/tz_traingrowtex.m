function [x,y] = tz_traingrowtex(img,cellbody,nucbody)
%TZ_TRAINGROWTEX Train texture growing model.
%   X = TZ_TRAINGROWTEX(IMG,CELLBOCY,NUCBODY)
%   
%   [X,Y] = TZ_TRAINGROWTEX(...)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function [x,y] = tz_traingrowtex(img,mask)
%OVERVIEW
%   
%PARAMETERS
%   img - 
%   cellbody - 
%   nucbody -
%RETURN
%   x - 
%   y - 
%DESCRIPTION
%   
%HISTORY
%   12-Apr-2005 Initial write TINGZ
%SEE ALSO
%   

imgsize=size(img);

if ~isempty(cellbody)    
    cellmask=tz_obj2img(cellbody,imgsize,{'2d','bn'});
    celledgeimg=bwperim(cellmask,8);
    celldistimg=bwdist(celledgeimg);
else
    cellmask=ones(imgsize);
    celldistimg=zeros(imgsize);
end

if ~isempty(nucbody)
    dnamask=tz_obj2img(nucbody,imgsize,{'2d','bn'});
    dnaedgeimg=bwperim(dnamask,8);
    dnadistimg=bwdist(dnaedgeimg);
else
    dnamask=zeros(imgsize); 
    dnadistimg=zeros(imgsize);
end

mask=cellmask-dnamask;

img=double(img);

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
zs=zeros(length(maskpos),1);
cellpos=celldistimg(maskpos);
dnapos=dnadistimg(maskpos);
x=[];%[1+zs,1+zs,cellpos,dnapos,zs];
y=[];%img(maskpos)==niter;

for k=niter-1:-1:1
    k
    prevsimg=(img>k);
    cursimg=(img==k);
    prevpkernel=pkernel;
    prevtkernel=tkernel;
    
    curmaskpos=maskpos;
    curmaskpos(prevsimg(maskpos)==1)=[];
    
    nones=nbimg(curmaskpos);
    cellpos=celldistimg(curmaskpos);
    dnapos=dnadistimg(curmaskpos);
    zs=zeros(length(nones),1);
    x=[x;[1+zs,niter-k+zs+1,cellpos,dnapos,nones]];
    y=[y;cursimg(curmaskpos)];
    
    nbimg=conv2(cursimg,nbkernel,'same')+nbimg;
end


