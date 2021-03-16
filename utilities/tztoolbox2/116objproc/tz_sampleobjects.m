function [selobjects,dists2] = ...
    tz_sampleobjects(combobjects,clstlabel,objnum,occuridx,dists)
%TZ_SAMPLEOBJECTS Sampling objects from clusters.
%   SELOBJECTS = 
%       TZ_SAMPLEOBJECTS(COMBOBJECTS,CLSTLABEL,OBJNUM,OCCURIDX)
%   returns a cell array of objects drawn from COMBOBJECTS randomly.
%   It is supposed that the objects have been sorted into some clusters,
%   which is described by the vector CLSTLABEL here. OBJNUM is the a
%   vector of object numbers from the clusters in OCCURIDX. SELOBJECTS
%   will be sorted by their sizes from large to small.
%
%   [SELOBJECTS,DISTS2] = 
%       TZ_SAMPLEOBJECTS(COMBOBJECTS,CLSTLABEL,OBJNUM,OCCURIDX,DISTS)
%   also specify the distances of the objects.
%   
%   See also

%   23-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 4
    error('4 or 5 arguments are required')
end

selectedIndices = [];
for i=1:length(occuridx)
    availableIndex = find(clstlabel==occuridx(i));
    allPermutation = randperm(length(availableIndex));
    selectedIndex = availableIndex(allPermutation(1:objnum(i)));
    selectedIndices=[selectedIndices; selectedIndex];
end
selectedObjects = combobjects(selectedIndices);
combobjsize = tz_calcobjpsize(selectedObjects);
[sizes,rankIndex] = sort(combobjsize);
selobjects = selectedObjects(rankIndex,'descend');
if exist('dists','var')
    dists2 = dists(rankIndex);
end