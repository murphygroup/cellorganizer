function coords=tz_calcobjcof_mcf(objects)
%TZ_CALCOBJCOF_MCF COFs of MCF objects.
%   COORDS = TZ_CALCOBJCOF_MCF(OBJECTS) returns a 2-column matrix of COFs
%   of the objects in the cell array OBJECTS. See TZ_CALCOBJCOF for more
%   details.

%   ??-???-2004 Initial write T. Zhao
%   03-NOV-2004 Modified T. Zhao
%       - add comments
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argument is required')
end

nclass = length(objects);

for i=1:nclass
    ncell=length(objects{i});
    classcoords={};
    for j=1:ncell
        nobj=length(objects{i}{j});
        cellcoords=[];
        for k=1:nobj
            cellcoords(k,:)=tz_calcobjcof(objects{i}{j}{k});
        end
        classcoords{j}=cellcoords;
    end
    coords{i}=classcoords;
end

