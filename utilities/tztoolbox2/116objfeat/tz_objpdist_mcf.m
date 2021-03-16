function allpdists=tz_objpdist_mcf(objects1,objects2,imgsize,distset,option)
%TZ_OBJPDIST_MCF Distances of object pixels for MCF.
%   OBJECTS2 = TZ_OBJPDIST(OBJECTS1,OBJECTS2,IMGSIZE,DISTSET,OPTION)
%   
%   [OBJECTS2,IMGSIZE,DISTSET,OPTION] = TZ_OBJPDIST(...)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function tz_objpdist_mcf(objects1,objects2,imgsize,distset,option)
%
%OVERVIEW:
%   calculate pixel distance to objects
%PARAMETERS:
%   objects1 - objects for calculaing distance
%   objects2 - reference objects
%   imgsize - image size
%   distset - 
%   option - 
%RETURN:
%   alldists
%DESCRIPTION:
%   
%HISTORY:
%   ??-OCT-2004 Initial write TINGZ
%   01-NOV-2004 Modified TINGZ
%       - add comments


nclass=length(objects1);

for i=1:nclass
    ncell=length(objects1{i});
    
    for j=1:ncell
        [i j]
        allpdists{i}{j}=tz_objpdist(objects1{i}{j},objects2{i}(j),imgsize,distset,option);
    end
end