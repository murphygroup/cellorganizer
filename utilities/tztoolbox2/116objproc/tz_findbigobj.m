function [bigobj,classidx,cellidx] = ...
    tz_findbigobj(combobjects,combclass,combcellidx,threshold)
%TZ_FINDBIGOBJ Extract objects with size above the threshold.
%   BIGOBJ = TZ_FINDBIGOBJ(COMBOBJECTS,COMBCLASS,COMBCELLIDX,THRESHOLD)
%   
%   [BIGOBJ,CLASSIDX,CELLIDX] = TZ_FINDBIGOBJ(...)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function bigobj=tz_findbigobj(combobjects,threshold)
%
%OVERVIEW:
%   find objects with size above the threshold
%PARAMETERS:
%   combobjects - combined objects
%   combclass - combined class label
%   combcellidx - combined cell idx
%   threshold - size threshold
%RETURN:
%   bigobj - objects, cell array
%   classidx - class label
%   cellidx - cell idx
%DESCRIPTION:
%
%HISTORY:
%   ??-???-???? Initial write TINGZ
%   04-NOV-2004 Modified TINGZ

nobj=length(combobjects);

k=1;

for i=1:nobj
    if size(combobjects{i},1)>=threshold
        bigobj{k}=combobjects{i};
        classidx(k)=combclass(i);
        cellidx(k)=combcellidx(i);
        k=k+1;
    end
end

    