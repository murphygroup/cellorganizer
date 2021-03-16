function A=tz_deformnucell(cellcode,ncomp,nangle)
%TZ_DEFORMNUCELL Estimate deformation between a cell its nucleus.
%   A = TZ_DEFORMNUCELL(CELLCODE,NCOMP,NANGLE) returns the transformation
%   matrix for the deformation from a cell to its nucleus. The cell is 
%   described by CELLCODE. NCOMP is the number of components for
%   trasformation. For example, if NCOMP is 2, then the nucleus will be
%   separated into two parts to match the cell, and A will be a 3x3x2
%   matrix. NANGLE specifies the angles parcipitating in transformation. If
%   it is an integer, the set of angles for transfomation will be at NANGLE
%   orientations with equal interval. If it is a vector, NANGLE itself will
%   be the set of angles.
%   
%   See also TZ_ESTTRANSFNUCELL

%   ??-???-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 3
    error('Exactly 3 arguments are required')
end

if length(nangle)==1
    [maxcelldist,initangle]=max(cellcode.nucdist);
    astep=360/nangle;
    selangle=linspace(1,360-astep,nangle);
    selangle=selangle-astep+initangle-1;
    selangle=round(selangle);
    selangle(selangle<=0)=selangle(selangle<=0)+360;
    selangle(selangle>360)=selangle(selangle>360)-360;
    psangle=ml_redistrnum(length(selangle),ncomp);
else
    selangle=nangle;
%     cangle=360/ncomp;
%     for i=1:ncomp
%         psangle(i)=sum(selangle<=i*cangle & selangle>(i-1)*cangle);
%     end
%     
%     if any(psangle<=2)
        psangle=tz_redistrnum(length(selangle),ncomp);
%     end
end


affpts=[];
nuccenter=mean(cellcode.nuchitpts,1);

tz_showpts_2d(ml_addrow(cellcode.cellhitpts,-nuccenter),'ln',1);
hold on

for i=1:ncomp
    avgangle=mean(selangle(1:psangle(i)))*pi/180;
    selcellpts=cellcode.cellhitpts(selangle(1:psangle(i)),:);
    selnucpts=cellcode.nuchitpts(selangle(1:psangle(i)),:);
    selnucpts=ml_addrow(selnucpts,-nuccenter);
    selnucpts=tz_rotate_2d(selnucpts,-avgangle);
    selcellpts=ml_addrow(selcellpts,-nuccenter);
    selcellpts=tz_rotate_2d(selcellpts,-avgangle);
    A(:,:,i)=tz_affinepts(selnucpts,selcellpts);
    tz_showpts_2d(ml_addrow(cellcode.nuchitpts(selangle(1:psangle(i)),:), ...
        -nuccenter),'ln',1);
    selangle(1:psangle(i))=[];
    
    % affpts=tz_transform_2d(tz_addrow(selnucpts,-nuccenter),A(:,:,1));  
    % affpts=[affpts;tz_transform_2d(tz_addrow(selnucpts,-nuccenter),A(:,:,2))];
end



% tz_showpts_2d(tz_addrow(cellcode.nuchitpts(selangle(5:8),:),-nuccenter),'ln',1);

% tz_showpts_2d(affpts,'cp',1);
% tz_showpts_2d(tz_addrow(cellcode.cellhitpts(selangle,:),-nuccenter),'ln',0);
hold off
