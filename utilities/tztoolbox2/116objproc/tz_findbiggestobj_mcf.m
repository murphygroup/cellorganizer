function [bigobj,classidx]=tz_findbiggestobj_mcf(objs)
%TZ_FINDBIGGESTOBJ_MCF Find biggest object in each cell of MCF.
%   BIGOBJ = TZ_FINDBIGGESTOBJ_MCF(OBJS) returns a cell array of objects,
%   which are the biggest objects of cells in the three-level cell array 
%   OBJS.
%   
%   [BIGOBJ,CLASSIDX] = TZ_FINDBIGGESTOBJ_MCF(...) also returns class
%   labels for BIGOBJ.

%   ??-???-???? Initial write T. Zhao
%   04-NOV-2004 Modified T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argument is required')
end

nclass=length(objs);

k=1;
bigobj={};
classidx=[];

for i=1:nclass
    ncell=length(objs{i});
    for j=1:ncell
        sizes=tz_objsizes(objs{i}{j});
        [ignore,p]=max(sizes);
        bigobj{k}=objs{i}{j}{p};
        classidx(k)=i;
        k=k+1;
    end
end
