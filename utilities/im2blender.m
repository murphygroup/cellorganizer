function [shiftvector, answer] = im2blender( img, savefile, downsample, patchsample, shiftvector_flag)
% IMG2BLENDER Exports a generated instance from CellOrganizer to a .obj mesh format
% that can be read by Blender.
%
%
% List Of Input Arguments     Descriptions
% -----------------------     ------------
% img                         a 3D image you wish to obtain the mesh for
% downsample                  factor by which you wish to downsample 
% savefile                    the path and file name you wish to save the generated file as
% patchsample                 the percentage of the verticies calculated that the user
%                             wants kept. Keeping more will result in a larger .obj file but have
%                             better resolution. Default value is 0.05
% shiftvector                 a 1x3 vector of the amount to shift the resulting mesh.
%
% This is used to place the mesh at the origin when called with syn2blender.
%
% List Of Outputs             Descriptions
% ---------------             ------------
% shiftvector                 a 1x3 vector of the amount to shift the resulting mesh.
% answer                      boolean flag that indicates if the test was successful

% Author: Devin Sullivan 
%
% Copyright (C) 2012-2016 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
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

% May 8, 2012 Devin Sullivan First downsample img to be on a reasonable scale 
% where 1 pixel = 1micron for mcell simulation. For standard cell organizer
% this is a downsampling of 5
% June 19, 2012 I. Cao-Berg Changed the method from ml_makeobjfile.m to img2blender.m,
% added license, documentation as well as checking of the input argument
% July 26 2012 D. Sullivan Create centered object files and reorient cell
% so that the bottom is on the horizontal plane
% October 29, 2012 D. Sullivan Changed the isosurface threshold from 1 to
% 0.5 to find the isosurface for the boolean images. 
% Nov 9, 2012 D. Sullivan Added a reducepatch call to significantly reduce
% the file size for the .obj files (e.g. 70MB to ~100KB). Added the
% reduction fraction to the list of inputs and an if statement to set a
% default if it's not defined.
% May 29, 2013 D. Sullivan added "shiftvector" parameter to maintain
% consistent shift (to origin) when running multiple objects in a single 

%icaoberg june 19, 2012
%5/29/13 D.Sullivan added checks for shiftvector
%step 0: check input arguments
answer = false;
if nargin == 3
    patchsample = 0.05;
    %    shiftvector = [0,0,0];
    % elseif nargin == 4
    %   shiftvector = [0,0,0];
elseif nargin < 3
    error('im2blender requires at least 3 input arguments');
end

if isempty( img )
    warning('Input argument image cannot be empty.')
    return
end

%5/29/13 D.Sullivan - added catch for empty patchsample
if isempty ( patchsample )
    patchsample = 0.05;
end

if ~isa( savefile, 'char' )
    warning('Input argument savefile must be a string')
    return
end

%need to resize this image
%must give the downsized ratio for each dimension
downsample3D = [downsample,downsample,downsample];
img = ml_downsize(img,downsample3D);

%need to pad the image so that the isosurface has no holes
img = padarray(img,[3,3,3]);


%7/26/12 D.Sullivan
D = smooth3(squeeze(img),'gaussian',[11,11,11]);
%D = smooth3(squeeze(img),'gaussian',[1,1,1]);

[x,y,z] = meshgrid(1:size(img,2),1:size(img,1),1:size(img,3));
x = x-mean(x(:));
y = y-mean(y(:));
z = z-mean(z(:));
v = D;

%make iso-surface (Mesh) of skin
%DPS 10/29/12 changed isosurface threshold for boolean images
%FV = isosurface(D,1);
% FV = isosurface(D,0.5);
%DPS 7/15/14 - changed isosurface threshold to be mean of image.
FV = isosurface(D,mean(D(:)));

% R. Arepally 6/7/13 checks if the shiftvector_flag is empty. If it is
% then calculate the mean of vertices. This becomes the shift vector for
% all of the objects that are created.
if isempty ( shiftvector_flag )
    shiftvector_flag = mean(FV.vertices);
end




%5/29/13 D.Sullivan changed from mean since mean of each object was slighly
%different and the meshes were no longer lined up.
% FV.vertices = FV.vertices-repmat(mean(FV.vertices),size(FV.vertices,1),1);
FV.vertices = FV.vertices-repmat(shiftvector_flag,size(FV.vertices,1),1);
%DPS 11/9/12 added reducepatch to significantly reduce the number of
%verticies calculated for the image
% FV = reducepatch(FV,0.05);
FV = reducepatch(FV,patchsample);
p = patch(FV);
set(p, 'FaceColor', 'red', 'EdgeColor', 'none');

%calculate iso-normals of the surface
N = isonormals(x,y,z,v,p);
L = sqrt(N(:,1).^2+N(:,2).^2+N(:,3).^2)+eps;
N(:,1) = N(:,1)./L; N(:,2)=N(:,2)./L; N(:,3)=N(:,3)./L;

% FV.faces = [FV.faces(:,3) FV.faces(:,2) FV.faces(:,1)];
material(1).type='newmtl';
material(1).data='skin';
material(2).type='Ka';
material(2).data=[0.8 0.4 0.4];
material(3).type='Kd';
material(3).data=[0.8 0.4 0.4];
material(4).type='Ks';
material(4).data=[1 1 1];
material(5).type='illum';
material(5).data=2;
material(6).type='Ns';
material(6).data=27;
clear OBJ
OBJ.vertices = FV.vertices;
OBJ.vertices_normal = N;
OBJ.material = material;
OBJ.objects(1).type='g';
OBJ.objects(1).data='skin';
OBJ.objects(2).type='usemtl';
OBJ.objects(2).data='skin';
OBJ.objects(3).type='f';
OBJ.objects(3).data.vertices=FV.faces;
OBJ.objects(3).data.normal=FV.faces;
write_wobj(OBJ,savefile);

answer = true;
shiftvector = shiftvector_flag;
