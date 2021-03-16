function cellpts=tz_expandnuc(nucbody,imgsize,angles,rdist)
%TZ_EXPANDNUC Extract points from nucleus edge at different angles.
%   CELLPTS = TZ_EXPANDNUC(NUCBODY,IMGSIZE,ANGLES,RDIST) returns an array
%   of points of cell boundary that is expanded from the nucleus NUCBODY,
%   which is also an array of points. There points correspond to the set of
%   angles specified by the vector ANGLES. RDIST is also a vector and it 
%   has the same size as ANGLES. Each element in RDIST is the ratio of the
%   distance from nucdist to celldist, where nucdist is the distance from
%   nuclear center to nuclear boundary and celldist is the distance from
%   nuclear center to cell boundary. IMGSIZE is the size of image for
%   implementation.
%   
%   See also TZ_EXPANDNUC2

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 4
    error('Exactly 4 arguments are required')
end

nucimg=tz_obj2img(nucbody,imgsize,{'2d','bn'});
nucedge=bwperim(nucimg,8);
nuccenter=round(mean(nucbody,1));
len=sqrt(sum(imgsize.^2));
nucdist=[];
nuchitpts=[];
if isempty(angles)
    for a=0:1:359
        pts=tz_getlinept2(nuccenter(1:2),a,len);
        ps=tz_imgptspixel(nucedge,pts);
        intc1=find(ps>0);
   
        if ~isempty(intc1)
            nuchitpts=[nuchitpts;pts(intc1(end),:)];
            nucdist=[nucdist,sqrt(sum((nuccenter-nuchitpts(end,:)).^2))];
            angles=[angles,a];
        else
            imshow(tz_setimglnpixel2(double(nucedge)+double(celledge), ...
                nuccenter(1:2),a+theta,len),[]);   
        end
    end
    
    [maxdist,theta]=max(nucdist);
    theta=theta-1;
    angles=mod(theta+(0:359),360);
end

dists=nucdist(angles+1);
celldist=dists./rdist;

for i=1:length(dists)
    pts=tz_getlinept2(nuccenter(1:2),angles(i),celldist(i),1);
    cellpts(i,:)=pts(end,:);
end




