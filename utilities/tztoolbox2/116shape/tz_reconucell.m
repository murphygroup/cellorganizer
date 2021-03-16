function tz_reconucell(cellcode,A,nangle)
%TZ_RECONUCELL Show transformation from nucleus to cell.
%   TZ_RECONUCELL(CELLCODE,A,NANGLE) plots the transformed nucleus and
%   original nucleus at the same time for CELLCODE. A is the transformation
%   matirx and NANGLE specifies the angles for transformation. It could be
%   a vector or an integer.
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function tz_reconucell(cellcode,A)
ncomp=size(A,3);
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
%     selangle=nangle;
%     cangle=360/ncomp;
%     for i=1:ncomp
%         psangle(i)=sum(selangle<=i*cangle & selangle>(i-1)*cangle);
%     end
%      if any(psangle<=1)
        psangle=tz_redistrnum(length(selangle),ncomp);
%     end
end


affpts=[];
nuccenter=mean(cellcode.nuchitpts,1);

for i=1:ncomp
    avgangle=mean(selangle(1:psangle(i)))*pi/180;
    selcellpts=cellcode.cellhitpts(selangle(1:psangle(i)),:);
    selnucpts=cellcode.nuchitpts(selangle(1:psangle(i)),:);
    selnucpts=ml_addrow(selnucpts,-nuccenter);
    selnucpts=tz_rotate_2d(selnucpts,-avgangle);
    selcellpts=ml_addrow(selcellpts,-nuccenter);
    selcellpts=tz_rotate_2d(selcellpts,-avgangle);
    temppts=tz_transform_2d(selnucpts,A(:,:,i));
    temppts=tz_rotate_2d(temppts,avgangle);   
    affpts=[affpts;temppts];
    selangle(1:psangle(i))=[];
end
% tz_showpts_2d(tz_addrow(cellcode.cellhitpts,-nuccenter),'ln',1);
hold on

tz_showpts_2d(ml_addrow(cellcode.nuchitpts,-nuccenter),'ln',1);
% tz_showpts_2d(tz_addrow(cellcode.nuchitpts(selangle(5:8),:),-nuccenter),'ln',1);

tz_showpts_2d(affpts,'cp',1);
% tz_showpts_2d(tz_addrow(cellcode.cellhitpts(selangle,:),-nuccenter),'ln',0);
hold off
