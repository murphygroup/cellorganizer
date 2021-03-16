function dists=tz_objpdist(objects1,objects2,imgsize,distset,option)
%TZ_OBJPDIST Distances from pixels of some objects to other objects.
%   OBJECTS2 = TZ_OBJPDIST(OBJECTS1,OBJECTS2,IMGSIZE,DISTSET,OPTION)
%   
%   [OBJECTS2,IMGSIZE,DISTSET,OPTION] = TZ_OBJPDIST(...)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function dists=tz_objpdist(objects1,objects2,imgsize,distset,option)
%OVERVIEW
%   Calculate distances between two objects
%PARAMETERS:
%   objects1 - the first object
%   objects2 - the second object, for building distance map
%   imgsize - size of the original image
%   distset - 'cof','edg'
%   option - signed or not
%RETURN:
%   dists - a cell array
%DESCRIPTION:
%   each pixel in objects1 is calculated
%
%HISTORY:
%   04-OCT-2004 Initial write TINGZ

objimg=tz_obj2img(objects2,imgsize)>0;
%%%%%%%%%%%%%%%%%%%%
% imshow(objimg,[]);
% drawnow
%%%%%%%%%%%%%%%

fobjimg=bwfill(objimg,'hole');
%%%%%%%%%%%%%%%%%
% imshow(fobjimg,[]);
% drawnow
%%%%%%%%%%%%%%%%%

dists={};

for i=1:length(distset)
    switch distset{i}
    case 'cof'
        distimg=zeros(imgsize);
        for j=1:length(objects2)
            cof=round(tz_calcobjcof(objects2{j}));
            distimg(cof(1),cof(2))=1;
        end
        distimg=bwdist(distimg);
        
    case 'edg'
        distimg=bwdist(objimg);
    end
    tmpdists=[];
    for j=1:length(objects1)
        ind=sub2ind(size(distimg),objects1{j}(:,1),objects1{j}(:,2));
        tmpdists=distimg(ind);
        if option(i)==1
            signs=fobjimg(ind);
            tmpdists(signs)=-tmpdists(signs);
        end
        dists{j}(:,i)=tmpdists;
    end
end
