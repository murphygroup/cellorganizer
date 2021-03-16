function [axlns,dists,borders,majorangles,centers] = ...
    tz_parsenucleus(combobjects,imgsize,isshow)
%TZ_PARSENUCLEUS Extract medial axis information from objects.
%   AXLNS = TZ_PARSENUCLEUS(COMBOBJECTS,IMGSIZE,ISMHOW)
%   
%   [AXLNS,DISTS,BORDERS,MAJORANGLES,CENTERS] = TZ_PARSENUCLEUS(...)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function [axlns,dists,borders]=tz_parsenucleus(combbojects,isshow)
%
%OVERVIEW:
%   extract medial axis information from objects
%PARAMETERS:
%   combobjects - object cell array
%   isshow - show the procedure or not
%RETURN:
%   axlns - medial axis
%   dists - width
%   borders - edge
%DESCRIPTION:
%   
%HISTORY:
%   ??-???-2004 Initial write TINGZ
%   05-NOV-2004 Modified TINGZ
%       - add comments

for i=1:length(combobjects)
    bwobj=combobjects{i};
    bwobj(:,3)=1;
    bwobjs{1}=bwobj;
    bwimg=tz_objs2img(bwobjs,imgsize,{'2d','bn'});
    bwimg2=tz_obj2img(bwobj,[],{'2d','bn'});
    
    imshow(bwimg,[]); 
    
    Trans=bwimg;
    [Theta,center] = tz_bwmajorangle(bwimg2);
    Theta=Theta*180./pi;
    majorangles(i)=Theta;
    centers(i,:)=center;
    
    Rot = imrotate(Trans, Theta, 'nearest', 'crop');
    
    if isshow
        figure(1)
    end
    imshow(Rot,[]);
    
    %imgedge=edge(opimg,'sobel');
    imgedge=bwperim(Rot);
    if isshow
        figure(2)
    end
    imshow(imgedge,[]);
    
    
    [imgaxis,axln,dist,border]=tz_imaxis(imgedge');
    [imgaxis2,axln2,dist2,border2]=tz_imaxis(imgedge);
    
    
    if size(axln2,1)>size(axln,1)
        imgaxis=imgaxis2;
        axlns{i}=axln2;
        dists{i}=dist2;
        borders{i}=border2;
    else
        axlns{i}=axln;
        dists{i}=dist;
        borders{i}=border;
    end
    
    imshow(imgaxis,[]);
    
    dist=dists{i};
    len=length(dist);
    if(mean(dist(1:floor(len/2)))>mean(dist(ceil(len/2):end)))
        dists{i}=fliplr(dists{i});
        axlns{i}=flipud(axlns{i});
        borders{i}=flipud(borders{i});
    end
    
    if isshow
        figure(3)
    end
    plot(dists{i});
    
%     title(num2str(i))
    i  
    drawnow
end

