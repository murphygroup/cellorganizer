function [feats, names] = tz_objskelfeats( objimg)
%TZ_OBJSKELFEATS Skeleton features for an object image.
%   FEATS = TZ_OBJSKELFEATS(OBJIMG)
%   
%   [FEATS,NAMES] = TZ_OBJSKELFEATS(...)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

% [FEATS, NAMES] = SKELFEATS( OBJIMG)
%
% Calculate skeleton features for the object OBJIMG.

objimg(:,sum(objimg,1)==0)=[];
objimg(sum(objimg,2)==0,:)=[];

objimg = double(objimg);
objbin = uint8(objimg>0);
%objskel = ml_mmthin_modify( objbin);
objskel = ml_mmthin(objbin);
skellen = length(find(objskel));
objsize = length(find(objimg));
skel_obj_area_ratio = skellen / objsize;
if sum(objskel(:))<3
    hullsize=skellen;
else
    skelhull = ml_imgconvhull( objskel);
    hullsize = length(find(skelhull));
    % if hull size comes out smaller than length of skeleton then it
    % is obviously wrong, therefore adjust
    if( hullsize < skellen) hullsize = skellen; end
end
skel_hull_area_ratio = skellen / hullsize;
skel_fluor = sum(objimg(find(objskel)));
obj_fluor = sum(objimg(:));
skel_obj_fluor_ratio = skel_fluor/obj_fluor;
branch_points = ml_find_branch_points(double(objskel));
no_of_branch_points = length(find(branch_points));
feats = [skellen skel_hull_area_ratio skel_obj_area_ratio ...
	 skel_obj_fluor_ratio no_of_branch_points/skellen];

names = {'obj_skel_len' ...
	 'obj_skel_hull_area_ratio' ...
	 'obj_skel_obj_area_ratio' ...
	 'obj_skel_obj_fluor_ratio' ...
	 'obj_skel_branch_per_len'};
