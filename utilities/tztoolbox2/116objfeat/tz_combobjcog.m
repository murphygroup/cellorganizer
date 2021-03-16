function centers=tz_combobjcog(combobjects)
%TZ_COMBOBJCOG Geometric center of objects
%   CENTER = TZ_COMBOBJCOG(COMBOBJECTS)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function tz_getobjcenter(combobjects)
%
%OVERVIEW:
%   calculate geometric center of each object of combined objects
%PARAMETERS:
%   combobjects - combined objects
%RETURN:
%   centers - object centers
%DESCRIPTION:
%
%HISTORY:
%   ??-???-???? Initial write TINGZ
%   04-NOV-2004 Modified TINGZ
%       - add comments
%       - remove flipping
%       - change funciton name tz_getobjcenter --> tz_combobjcog

for i=1:length(combobjects)
    centers(i,:)=mean(combobjects{i}(:,1:2),1);
end

%centers=fliplr(centers);