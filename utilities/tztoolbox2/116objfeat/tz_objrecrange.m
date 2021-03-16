function range=tz_objrecrange(obj)
%TZ_OBJRECRANGE Obsolete. See ML_OBJRECRANGE.
%   RANGE = TZ_OBJRECRANGE(OBJ)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function range=tz_objrecrange(obj)
%
%OVERVIEW:
%   find the size of the rectangle containing the object
%PARAMETERS:
%   obj - object 3xn
%RETURN:
%   range - [width,height]
%DESCRIPTION:
%
%HISTORY:
%   27-JUN-2004 Initial write TINGZ
%   04-NOV-2004 Modified TINGZ
%       - change funcion name tz_findrange --> tz_objrecrange

error(tz_genmsg('of','TZ_OBJRECRANGE','ML_OBJRECRANGE'));

% range(1)=max(obj(:,1))-min(obj(:,1))+1;
% range(2)=max(obj(:,2))-min(obj(:,2))+1;
range=max(obj,[],1)-min(obj,[],1)+1;