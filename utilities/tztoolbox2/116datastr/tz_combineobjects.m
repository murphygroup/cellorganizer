function [combobjects,combcellidx,combclass,combobjidx] ...
    = tz_combineobjects(objects,all_features)
%TZ_COMBINEOBJECTS Combine mcf objects.
%   COMBOBJECTS = TZ_COMBINEOBJECTS(OBJECTS,ALL_FEATURES) returns the 
%   combined objects, which is a one level cell array. OBJECTS has three
%   levels.
%   
%   COMBOBJECTS = TZ_COMBINEOBJECTS(OBJECTS,ALL_FEATURES) will also help
%   checking if OBJECTS and ALL_FEATURES have matched size.  ALL_FEATURES 
%   is the two-level cell array of object feature matrices.
%
%   [COMBOBJECTS,COMBCELLIDX,COMBCLASS,COMBOBJ] = TZ_COMBINEOBJECTS(...)
%   has not been available.

%   ??-???-???? Initial write T. Zhao
%   03-NOV-2004 Modified T. Zhao
%       - add comments (objects,all_features) --> (all_features,objects)
%       - change parameters 
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('1 or 2 arguments are required')
end

%warning('updated on 03-NOV-2004')

if exist('all_features','var')
    nclass=length(all_features);
else
    nclass = length(objects);
end

combobjects={};
combcellidx = [];
combclass = [];
combobjidx = [];

nobj=1;
for i=1:nclass
    if exist('all_features','var')    
        ncell=length(all_features{i});
    else
        ncell = length(objects{i});
    end

    for j=1:ncell
        combobjects = {combobjects{:},objects{i}{j}{:}};
        combclass = [combclass; zeros(length(objects{i}{j}),1)+i];
        combcellidx = [combcellidx; zeros(length(objects{i}{j}),1)+j];
        combobjidx =  [combobjidx; (1:length(objects{i}{j}))'];
        
%         for k=1:length(objects{i}{j})
%             combobjects{nobj}=objects{i}{j}{k};
%             nobj=nobj+1;
%         end
        if exist('all_features','var')
            if(length(objects{i}{j})~=size(all_features{i}{j},1))
                warning('cell unmatched');
            end
        end
    end
end
