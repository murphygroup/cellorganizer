function sizes=tz_objsizes(combobjs)
%TZ_OBJSIZES Number of pixels of each object in a cell array.
%   SIZES = TZ_OBJSIZES(COMBOBJS)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function sizes=tz_objssizes(combobjs)
%
%OVERVIEW:
%   get size of each object
%PARAMETERS:
%   combobjs - object cell array
%RETURN:
%   sizes - size vector
%DESCRIPTION:
%
%HISTORY:
%   ??-???-2004 Initial write TINGZ
%   04-NOV-2004 Modified TINGZ
%       - add comments
%       - change function name tz_getobjsize --> tz_objsizes

for i=1:length(combobjs)
    sizes(i)=size(combobjs{i},1);
end