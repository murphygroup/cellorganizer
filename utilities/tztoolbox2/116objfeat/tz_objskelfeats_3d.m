function  [feats, names] = tz_objskelfeats_3d(objects)
%TZ_OBJSKELFEATS_3D Skeleton features for an 3D object.
%   FEATS = TZ_OBJSKELFEATS_3D(OBJECTS)
%   
%   [FEATS,NAMES] = TZ_OBJSKELFEATS_3D(...)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function tz_objskelfeats_3d(objects)
%
%OVERVIEW:
%   calculate skeleton features for objects
%PARAMETERS:
%   objects - 3d objects
%RETURN
%   feats - skeleton features
%   names - name of the features
%DESCRIPTION:
%   the features are calculated on the slice with most above-threshold pixels
%
%HISTORY:
%   27-NOV-2004 Initial write TINGZ

for i=1:length(objects)
    object=objects{i}.voxels';
    z1=min(object(:,3));
    z2=max(object(:,3));
    sind=z1;
    maxsp=0;
    for j=z1:z2
        sp=sum(object(:,3)==j);
        if sp>maxsp
            maxsp=sp;
            sind=j;
        end
    end
    objimg=tz_obj2img(object(object(:,3)==sind,1:2));
    [feats(i,:),names]=tz_objskelfeats(objimg);
end