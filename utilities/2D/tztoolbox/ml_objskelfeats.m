function [feats, names] = ml_objskelfeats( objimg)

% [FEATS, NAMES] = ML_OBJSKELFEATS( OBJIMG)
% 
% Calculate skeleton features for the object OBJIMG.

% Copyright (C) 2006  Murphy Lab
% Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.
%
% For additional information visit http://murphylab.web.cmu.edu or
% send email to murphy@cmu.edu

% Written by Meel Velliste

%tz- 23-Jun-2006
% objimg = double(objimg);
%tz--

%tz+ 23-Jun-2006
objimg = ml_imcropbg(double(objimg));
%tz++

objbin = uint8(objimg>0);
% objskel = mmthin( objbin);
objskel = ml_mmthin( objbin);
% rewrite the mmthin function by yenixsa and Sam, 09/14/2004
skellen = length(find(objskel));
objsize = length(find(objimg));
skel_obj_area_ratio = skellen / objsize;
skelhull = ml_imgconvhull( objskel);
hullsize = length(find(skelhull));
% if hull size comes out smaller than length of skeleton then it
% is obviously wrong, therefore adjust
if( hullsize < skellen) hullsize = skellen; end
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
    
