function [combobjects,combcellidx,combclass]=tz_mcfobj2comb(objects)
%TZ_MCF2COMBOBJ Convert mcf objects to combined objects.
%   COMBOBJECTS = TZ_MCF2COMBOBJ(OBJECTS) returns the one-level cell
%   array of objects, which is the combination of the three-level cell
%   array of objects.
%
%   [COMBOBJECTS,COMBCELLIDX,COMBCLASS] = TZ_MCF2COMBOBJ(...) also
%   returns cell indices and class labels.

%   ??-???-2004 Initial write T. Zhao
%   03-NOV-2004 Modified T. Zhao
%       - add comments
%       - change function name tz_combineobjects_mcf --> tz_mcfobj2comb
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argument is required')
end

nclass=length(objects);
combobjects={};
combclass=[];
combcellidx=[];
nobj=1;
for i=1:nclass
    ncell=length(objects{i}); 
    for j=1:ncell
        combobjects={combobjects{1:end},objects{i}{j}{1:end}};
        combclass=[combclass;zeros(length(objects{i}{j}),1)+i];
        combcellidx=[combcellidx;zeros(length(objects{i}{j}),1)+j];
        nobj=nobj+1;
    end
end