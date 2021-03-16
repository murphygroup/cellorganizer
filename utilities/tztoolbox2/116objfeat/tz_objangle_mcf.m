function cas=tz_objangle_mcf(objcofs,cofs,ras)
%TZ_OBJANGLE_MCF Calculate angles of objects in cells (mcf)
%   CAS = TZ_OBJANGLE_MCF(OBJCOFS,COFS,RAS)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function cas=tz_objangle_mcf(objcofs,cofs,ras)
%
%OVERVIEW:
%   Calculate angles of objects in cells (mcf)
%PARAMETERS:
%   objcofs - cof of objects
%   cofs - cof of nucleus
%   ras - orientation of major axis of the cell
%RETURN:
%   cas - normalzied angles
%DESCRIPTION:
%
%HISTORY
%   ??-???-???? Initial write TINGZ
%   31-OCT-2004 Modified TINGZ
%       - add comments
%       - change function name tz_cellangle_mcf-->tz_objangle_mcf
%       
nclass = length(objcofs);

for i=1:nclass
    ncell=length(objcofs{i});
    ca={};
    for j=1:ncell
        ca{j}=tz_cellangle(objcofs{i}{j},cofs{i}(j,:),ras{i}(j,:));
    end
    cas{i}=ca;
end