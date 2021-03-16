function [beta,llk] = ...
    tz_trainlogitgrowtex(img,option,cellbody,nucbody,gsdir)
%TZ_TRAINLOGITGROWTEX Estimate logistic texture growing model.
%   BETA = TZ_TRAINLOGITGROWTEX(IMG,OPTION,CELLBODY,NUCBODY)
%   
%   BETA = TZ_TRAINLOGITGROWTEX(IMG,OPTION,CELLBODY,NUCBODY,GSDIR)
%   
%   [BETA,LK] = TZ_TRAINLOGITGROWTEX(...)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function beta = tz_trainlogitgrowtex(img,option,cellbody,nucbody,gsdir)
%OVERVIEW
%   ananlyze an image by growing texture on logistic regression
%PARAMETERS
%   img - input image
%   option - logistic data transformation
%   cellbody - cell body 
%   nucbody - dna body
%   gsdir - directiory for growing stages
%RETURN
%   beta - logistic model parameters
%DESCRIPTION
%   offset,time,celldist,nucdist,neighbor
%HISTORY
%   27-Apr-2005 Initial write TINGZ
%SEE ALSO
%   tz_traingrowtex

if ~exist('gsdir','var')
    gsdir=[];
end

imgsize=size(img);

if ~isempty(cellbody)    
    cellmask=tz_obj2img(cellbody,imgsize,{'2d','bn'});
    celledgeimg=bwperim(cellmask,8);
    celldistimg=bwdist(celledgeimg);
else
    cellmask=ones(imgsize);
end

if ~isempty(nucbody)
    dnamask=tz_obj2img(nucbody,imgsize,{'2d','bn'});
    dnaedgeimg=bwperim(dnamask,8);
    dnadistimg=bwdist(dnaedgeimg);
else
    dnamask=zeros(imgsize); 
end

mask=cellmask-dnamask;

img=double(img);

maskpos=find(mask==1);
img(find(mask==0))=0;
niter=max(img(:));

imgcode(:,1)=niter-img(maskpos)+1;
if ~isempty(cellbody)
    imgcode(:,2)=celldistimg(maskpos);
end

if ~isempty(nucbody)
    imgcode(:,3)=dnadistimg(maskpos);
end

maxcelldist=max(imgcode(:,2));

%8 directions
tmpimg1=[img(:,end),img(:,1:end-1)];
imgcode(:,4)=tmpimg1(maskpos);

tmpimg2=[img(:,2:end),img(:,1)];
imgcode(:,5)=tmpimg2(maskpos);

tmpimg3=[img(end,:);img(1:end-1,:)];
imgcode(:,6)=tmpimg3(maskpos);

tmpimg4=[img(2:end,:);img(1,:)];
imgcode(:,7)=tmpimg4(maskpos);

tmpimg=[tmpimg1(end,:);tmpimg1(1:end-1,:)];
imgcode(:,8)=tmpimg(maskpos);

tmpimg=[tmpimg2(end,:);tmpimg2(1:end-1,:)];
imgcode(:,9)=tmpimg(maskpos);

tmpimg=[tmpimg1(2:end,:);tmpimg1(1,:)];
imgcode(:,10)=tmpimg(maskpos);

tmpimg=[tmpimg2(2:end,:);tmpimg2(1,:)];
imgcode(:,11)=tmpimg(maskpos);
    
imgcode(:,4:11)=niter-imgcode(:,4:11)+1;

maxiter=100;
mine=1e-5;

%make zeros with the same size as training data
x=tz_mapdata(zeros(1,5),option);

if isempty(cellbody)
    x(:,3)=[];
    if isempty(nucbody)
        x(:,3)=[];
    end
else
    if isempty(nucbody)
        x(:,4)=[];
    end
end


beta=zeros(size(x,2),1);
countis=0;
llk=0;
for i=1:maxiter
    oldbeta=beta;
    delta=zeros(size(beta));
    H=zeros(length(beta),length(beta));
    oldllk=llk;
    llk=0;
%     tmpx=[];
%     tmpy=[];
%     tmpsx=[];
%     
    for j=2:niter
        curmaskpos=maskpos(imgcode(:,1)>=j);
        curimgcode=imgcode(imgcode(:,1)>=j,:);
        if ~isempty(curimgcode)
            y=curimgcode(:,1)==j;
            %         x=curimgcode(:,2:end);
            
            %x2: [cell distance, nuc distance]
            x2=curimgcode(:,2:3);
            
            %x3: neighbor weights
            if ~isempty(gsdir)
                load([gsdir '/' 'growstage' num2str(j-1) '.mat']);
                x3=growstageimg(curmaskpos);
            else
                x3=sum(curimgcode(:,4:end)<j,2);
            end
            %             px=[x2,x3];
            
            %x: [offset, sequence index, celldist, nucdist, neighbor weights]
            if strcmp(option,'gt11b')
                x2=x2/maxcelldist;
            end
            x=tz_mapdata([ones(length(y),1),ones(length(y),1)*j,x2,x3],option);
            
            if isempty(cellbody)
                x(:,3)=[];
                if isempty(nucbody)
                    x(:,3)=[];
                end
            else
                if isempty(nucbody)
                    x(:,4)=[];
                end
            end
            
%             px=x(:,3:end);
            sx=1./(1+exp(-(x*beta)));
%             delta(1)=delta(1)+sum(y-sx);
%             delta(2)=delta(2)+sum(y-sx)*j;
            dsx=y-sx;
%             for k=1:size(x,2)
%                 delta(k)=delta(k)+sum((dsx).*x(:,k));
%             end
            delta=delta+x'*dsx;
            csx=sx.*(1-sx);
%             H(1,:)=H(1,:)+csx'*x;
%             H(2,:)=H(2,:)+j*csx'*x;
            for k=1:size(x,2)
                H(k,:)=H(k,:)+(x(:,k).*csx)'*x;
            end
            
            llk=llk-sum((x(y==0,:)*beta))+sum(log(sx));
            countis=countis+sum(x(y==0,3));
%             countis=max([countis;y(x(:,end)==0)]);%max([countis;x(y==1,end)]);
            
%             tmpx=[tmpx;x];
%             tmpy=[tmpy;y];
%             tmpsx=[tmpsx;sx];
        end
    end
    
    if abs(oldllk-llk)<mine | isnan(llk) | isinf(llk)
        break
    end
    llk
    iH=inv(H);
    beta=beta+iH*delta
    sum(countis);
    if mean(abs(beta-oldbeta))<mine
        break
    end
   
end


