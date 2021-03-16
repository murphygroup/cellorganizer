function combobj2 = tz_combineobjcoords(all_coords,all_features)
%TZ_COMBINEOBJCOORDS Check and combine all object coodinates.
%   COMBOBJ2 = TZ_COMBINEOBJCOORDS(ALL_COORDS,ALL_FEATURES) returns the
%   the combined object coordinates in the MCF cell array ALL_COORDS.
%   It will also check if ALL_COORDS and ALL_FEATURES have the same size.

%   ??-???-???? Initial write T. Zhao
%   03-NOV-2004 Modified T. Zhao
%       - add comments (all_features,all_coord) --> all_coords,all_features
%       - change parameters 
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

warning('updated on 03-NOV-2004')

nclass=length(all_features);
combobj2=[];
for i=1:nclass
    ncell=length(all_features{i});
    for j=1:ncell
        combobj2=[combobj2;all_coords{i}{j}];
        if(size(all_coords{i}{j},1)~=size(all_features{i}{j},1))
            warning('cell unmatched');
        end
    end
end