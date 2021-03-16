function objects=tz_combobj2mcf(combobjects,combclass,combcellidx)
%TZ_COMBOBJ2COF Convert combined objects to MCF.
%   OBJECTS = TZ_COMBOBJ2COF(COMBOBJECTS,COMBCLASS,COMBCELLIDX) returns
%   the three-level cell array of objects. COMBOBJCTS is the one-level
%   cell array of the objects. COMBCLASS is the vector of class labels
%   and combcellidx is the vector of cell indices.

%   ??-???-???? Initial write T. Zhao
%   31-OCT-2004 Modified T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University


if nargin < 3
    error('Exactly 3 arguments are required')
end

classes=tz_findclass(combclass);

for i=1:length(classes)
    cellidx=combcellidx(combclass==classes{i}(1));
    ncell=max(cellidx);
    cellobjs=combobjects(combclass==classes{i}(1));
    for j=1:ncell
        objects{i}{j}=cellobjs(cellidx==j);
    end   
end