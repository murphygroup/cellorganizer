function newobjects = tz_rotateobjs(combobjects,theta,orgcenter,newcenter)
%TZ_ROTATEOBJS Rotate and translate objects. (Under construction)
%   NEWOBJECTS = TZ_ROTATEOBJS(COMBOBJECTS,THETA,ORGCENTER,NEWCENTER)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function newobjects = tz_rotateobjs(combobjects,theta,orgcenter,newcenter)
%
%OVERVIEW:
%   rotate and translate objects
%PARAMETERS:
%   combobjects - objects
%   theta - rotate angle
%   orgcenter - rotate center
%   newcenter - translated center
%RETURN:
%   newobject - new objects
%DESCRIPTION
%
%HISTORY:
%   ??-???-???? Initial write TINGZ
%   05-NOV-2004 Modified TINGZ
%       - add comments
%       - change function name tz_rotateobj --> tz_rotateobjs

for i=1:length(combobjects)
    newobjects{i}=[combobjects{i}(:,1)-orgcenter(1), ...
            combobjects{i}(:,2)-orgcenter(2)];
    newobjects{i}=tz_rotate_2d(newobjects{i},theta);
    newobjects{i}=[newobjects{i}(:,1)+newcenter(1), ...
            newobjects{i}(:,2)+newcenter(2),combobjects{i}(:,3)];
end

