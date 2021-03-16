function cas=tz_pangle_mcf(objects,cofs,ras)
%TZ_PANGLE_MCF Cell related angle for each pixel in objects.
%   CAS = TZ_PANGLE_MCF(OBJECTS,COFS,RAS)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function cas=tz_cellangle_mcf(objcofs,cofs,ras)
%
%OVERVIEW:
%   calculate cell related angle for each pixel in objects
%PARAMTERS:
%   objects - objects, cell array
%   cofs - cof of the nucleus
%   ras - angle of cell major axis
%RETURN:
%   cas - normalized angle
%
%HISTORY:
%   ??-OCT-2004 Initial write TINGZ
%   01-NOV-2004 Modified TINGZ

nclass = length(objects);

for i=1:nclass
    ncell=length(objects{i});
    
    for j=1:ncell
        [i j]
        nobj=length(objects{i}{j});
        ca={};
        for k=1:nobj
            ca{k}=tz_cellangle(objects{i}{j}{k}(:,1:2),cofs{i}(j,:),ras{i}(j,:));
        end
        cas{i}{j}=ca;
    end
   
end