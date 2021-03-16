function feats = ml_objfeat3d(objects,dnaobjs,imgsize,res)
%ML_OBJFEAT3D Calculate 3D object level features
%   FEATS = ML_OBJFEAT3D(OBJECTS,DNAOBJS,RES) returns a [feature matrix]
%   that are the features of the cell array of [3D object]s OBJECTs from 
%   the same cell. DNAOBJS is the [3D object]s of nucleus. RES is the
%   resolution, which is a 1x3 vector meaning [x y z].
%   
%   The calculated features are:
%       3D-SOF1.1: Number of voxels in object
%       3D-SOF1.2: 3D Euclidean distance between obj COF and DNA COF
%       3D-SOF1.3: Fraction of object voxels overlapping with DNA
%       3D-SOF1.4: A measure of eccentricity of the object
%       3D-SOF1.5: Number of holes in the object
%       3D-SOF1.6: A measure of roundness of the object
%       3D-SOF1.7: The length of the object's skeleton by homotopic thinning
%       3D-SOF1.8: The ratio of skeleton length to the area of the convex hull of the skeleton
%       3D-SOF1.9: The fraction of object pixel contained within the skeleton
%       3D-SOF1.10: The fraction of object fluorescence contained within the skeleton
%       3D-SOF1.11: The ratio of the number of branch points in skeleton to length of skeleton
%       3D-SOF1.12: Horizontal 2D Euclidean distance between object COF and DNA COF
%       3D-SOF1.13: Vertical 1D signed distance between object COF and DNA COF  
%   
%   See also

%   05-Dec-2006 Initial write T. Zhao
%   Copyright (c) 2006 Murphy Lab
%   Carnegie Mellon University
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation; either version 2 of the License,
%   or (at your option) any later version.
%   
%   This program is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   General Public License for more details.
%   
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
%   02110-1301, USA.
%   
%   For additional information visit http://murphylab.web.cmu.edu or
%   send email to murphy@cmu.edu


if nargin < 3
    error('Exactly 3 arguments are required');
end

dnaobj = [];
for i=1:length(dnaobjs)
    dnaobj = [dnaobj;dnaobjs{i}.voxels]; 
end

dnaimg=ml_obj2img(dnaobj,imgsize,{'3d','bn'});
dnacof = sum([dnaobj(:,1).*dnaobj(:,4),dnaobj(:,2).*dnaobj(:,4), ...
    dnaobj(:,3).*dnaobj(:,4)],1)/sum(dnaobj(:,4));
dnaind=sub2ind(imgsize,dnaobj(:,1),dnaobj(:,2),dnaobj(:,3));

for i=1:length(objects)
    disp([num2str(i) '/' num2str(length(objects))]);
    object=objects{i}.voxels;
    objind=sub2ind(imgsize,object(:,1),object(:,2),object(:,3));
    
    %3D-SOF1.1: Number of voxels in object
    feats(i,1)=size(object,1);
    
    %3D-SOF1.2: 3D Euclidean distance between obj COF and DNA COF
    %object cof
    objcof(i,:)=sum([object(:,1).*object(:,4),object(:,2).*object(:,4), ...
        object(:,3).*object(:,4)],1)/sum(object(:,4));
    %calculate later
    feats(i,2)=0; 
    
    objimg = ml_obj2img(object,[],{'3d','og'});
    
    %3D-SOF1.3: Fraction of object voxels overlapping with DNA
    objind(objind>length(dnaimg(:))) = [];
    feats(i,3)=sum(dnaimg(objind));
    
    objproj = sum(objimg,3);
    tmpobj = ml_findmainobj_bw(objproj>0);
    tmpimg=ml_obj2img(tmpobj,[]);
    
    %3D-SOF1.4: A measure of eccentricity of the object
    obj_feature = regionprops(tmpimg,'Eccentricity');
    feats(i,4) = obj_feature.Eccentricity;
    
    %3D-SOF1.5: Number of holes in the object
    feats(i,5) = objects{i}.n_holes;
    
    %3D-SOF1.6: A measure of roundness of the object
    objimageperim = bwarea(bwperim(tmpimg>0)) ;
    feats(i,6) = (objimageperim^2)/(4*pi*feats(i,1)) ;
    
    %3D-SOF1.7, 3D-SOF1.8, 3D-SOF1.9, 3D-SOF1.10, 3D-SOF1.11: Skeleton
    feats(i,7:11)=tz_objskelfeats(tmpimg);
    
    %3D-SOF1.12: Horizontal 2D Euclidean distance
    feats(i,12)=sqrt(sum(((objcof(i,1:2)-dnacof(1:2)).*res(1:2)).^2));  
    
    %3D-SOF2.13: Vertical 1D signed distance
    feats(i,13)=(objcof(i,3)-dnacof(3))*res(3);
    
end

feats(:,2)=ml_eucdist(dnacof',objcof',res)';
