function cas=tz_objangle_mcf(objcofs,cofs,ras)
%TZ_OBJANGLE_MDF Obsolete.

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

error(tz_genmsg('of','tz_objangle_mcf','tz_objangle_mcf'));

nclass = length(objcofs);

for i=1:nclass
    ncell=length(objcofs{i});
    ca={};
    for j=1:ncell
        ca{j}=tz_cellangle(objcofs{i}{j},cofs{i}(j,:),ras{i}(j,:));
    end
    cas{i}=ca;
end