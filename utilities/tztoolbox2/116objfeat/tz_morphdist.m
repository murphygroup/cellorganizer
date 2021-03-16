function [normdists,nucdists,celdists] = ...
    tz_morphdist(lcofs,celedge,nucedge,option)
%TZ_MORPHDIST Calculate normalized distance of object cof.
%   NORMDISTS = TZ_MORPHDIST(LCOFS,CELLEDGE,NUCEDGE,OPTION)
%   
%   [NORMDISTS,NUCDISTS,CELLDISTS] = TZ_MORPHDIST(...)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function normdist=tz_morphdist(lcof,celedge,nucedge,option)
%
%OVERVIEW:
%   Calculate normalized distance of object cof
%PARAMETERS:
%   lcof - matrix for labeled cof, nx3
%   celedge - cell edge image
%   nucedge - dna edge image
%   option - option for different scale schemes
%       'cn' - cell to nuclear
%       'mc' - cell to nuclear, nuclear to center
%RETURN:
%   normdists - normalized distance to DNA edge
%DESCRIPTION:
%   dnadist/(celldist+dnadist)
%
%HISTORY:
%   28-SEP-2004 Initial write TINGZ
%   29-SEP-2004 Modified TINGZ
%       - Add option to deal with different scale schemes


nobj=size(lcofs,1);

subplot(2,1,1)
imshow(nucedge,[]);
nucedge2=double(nucedge);
nucedge=tz_mainobjimg_bw(nucedge);
[r,c]=find(nucedge==1);
nucobj=[r,c];
subplot(2,1,2)
imshow(nucedge,[]);
drawnow
celdistimg=bwdist(celedge);
nucdistimg=bwdist(nucedge);

if strcmp(option,'mc')
    nuccof=mean(nucobj);
end

for i=1:nobj
    cof=round(lcofs(i,1:2));
    celdist=celdistimg(cof(1),cof(2));
    nucdist=nucdistimg(cof(1),cof(2));
    if lcofs(i,3)==1
        nucdist=-nucdist;
    end
    
    celdists(i)=celdist;
    switch option
    case 'cn'
        nucdists(i)=nucdist;
        if celdist+nucdist==0
            nucdist
        end
        normdists(i)=nucdist/(celdist+nucdist);
    case 'mc'
        nuccofdist=sqrt(sum((lcofs(i,1:2)-nuccof).^2));
        nucdists(i,1)=nucdist;
        nucdists(i,2)=nuccofdist;
        if nucdist>0
            normdists(i)=nucdist/(celdist+nucdist);
        else
            normdists(i)=nucdist/(nuccofdist+abs(nucdist));    
        end
    end
end