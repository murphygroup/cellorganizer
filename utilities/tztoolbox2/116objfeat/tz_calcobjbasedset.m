function feature=tz_calcobjbasedset(combobj,featset)
%TZ_CALCOBJBASEDSET Calculate object-based features.
%   FEATURE = TZ_CALCOBJBASEDSET(COMBOBJ,FEATSET) returns a row vector of
%   cell features. The available feature sets are:
%   (OBJCOM: object feature matrix)
%       'objnum' - object number (only depends on the number of rows of 
%           the feature matrix)
%       'featsum' - sum of the object features
%       'featmean' - mean of the object features
%   (OBJCOM: a vector of fluorescence fraction)
%       'fluofrac' - fluorescence fraction
%   (OBJCOM: a vector of fluorescence)
%       'totalfluo' - total fluorescence
%   (OBJCOM: a vector of object sizes)
%       'objsize' - object size
%   (OBJCOM: a vector of object distances)
%       'objdist' - mean and variance of object distance
%   (OBJCOM: a 2-column matrix of object center coordinates)
%       'mst' - minimal spanning tree features (for distance only)
%       'mst2','fmst' - the other set of MST features. See 
%           TZ_CALCTREEFEAT for more details.
%   (OBJCOM: a vector MST features)
%       'readymst' - just return COMBOBJ
%   (OBJCOM: a cell array of objects)
%       'pixelnum' - total number of pixels in all objects
%       'fluosum' - total fluorescence in all objects
%       
%   Notice:
%       Their meanings are really dependent on the input COMOBJ,
%       which is explained in the paranthesis.
%       This function is not well designed. 
%
%   See also

%   06-JUN-2004  Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

switch(featset)
case 'featsum'
    feature=sum(combobj,1);
case 'featmean'
    feature=mean(combobj,1);
case 'objnum'
    feature=size(combobj,1);
case {'fluofrac','totalfluo','objsize'}
    feature=sum(combobj);
case 'objdist'
    if isempty(combobj)
        feature=[0 0];
    else
        feature=[mean(combobj) var(combobj)];
    end
case 'mst'
    if size(combobj,1)<=1
        feature=zeros(1,2);
    else
        E=tz_mstree_l(combobj);
        feature=[mean(E(3,:)),var(E(3,:))];
    end
case {'mst2','fmst'}
    if size(combobj,1)<=3
        E=tz_mstree(combobj,'eu',1);
    else
        E=tz_mstree_l(combobj);
    end
    feature=tz_calctreefeat(E',combobj);
case 'readymst'
    feature=combobj;
case 'pixelnum'
    feature=0;
    for i=1:length(combobj)
        feature=feature+size(combobj{i},1);
    end
case 'fluosum'
    feature=0;
    for i=1:length(combobj)
        feature=feature+sum(combobj{i}(:,3));
    end
otherwise
    feature=[];
end