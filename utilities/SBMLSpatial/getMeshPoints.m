function FV = getMeshPoints(image, savepath, downsampling, options)
%GETMESHPOINTS This function gets mesh points given an image
%
%Inputs:
%image = a 2D or 3D image
%savepath = path to save image after preprocessing (D.mat), can be empty
%downsampling = downsampling fraction (1 means no downsampling, 1/5 means 1/5 the size)
%
%Outputs:
%FV = a struct containing the vertices etc for the mesh

% Author: Devin Sullivan September 19,2013
% Author: Taraz Buck January 29, 2019
%
% Copyright (C) 2012-2019 Murphy Lab
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

warning('off', 'CellOrganizer:getMeshPoints');

if nargin<2 || isempty(savepath)
    savepath = pwd;
end
if nargin<3
    downsampling = 1;
end
if length(downsampling)==1
    downsampling = downsampling * [1,1,1];
end

%Parameters
%%%
padAmount = 3;%This creates padding to close the sides, top and bottom of the mesh
%ideally this would be resolution dependent
blurAmount = 11;%This smooths the image slightly to create a more smoothed mesh and reduce artifacts from pixelization
%This value was empirically chosen for the HeLa data, but may be different
%for different membrane properties, resolutions etc
%should not hard code this.
patchsample = 0.05;%This is the amount of sub-sampling to be done.
warning('CellOrganizer:getMeshPoints', 'Shouldn''t the threshold default to 0.5 if the image is binary?')
isopatch = 0.9;
%%%

disp('Entering mesh code.')
warning('CellOrganizer:getMeshPoints', 'Converting image to single due to memory issues');
image = single(image);
%First trim the image to ensure we will not mesh things we don't need to.
%image2 = image;
xvals = find(sum(sum(image,3),1));
yvals = find(sum(sum(image,3),2));
zvals = find(sum(sum(image,1),2));
disp('Trimming image')
image = image(min(yvals):max(yvals),min(xvals):max(xvals),min(zvals):max(zvals));
%Keep track of how much is cut off by this

disp('Padding image')
%need to pad the image so that the isosurface has no holes
image = padarray(image,[padAmount,padAmount,padAmount]);

%Smooth the image a little, should probably not hard code 11.
disp('Smoothing image')
image_size = size(image);
smoothing_time = tic;
% This uses smooth3's default sigma of 0.65
% D = smooth3(squeeze(image),'gaussian',[blurAmount,blurAmount,blurAmount]);
% Runs out of memory on lanec1 model1
% D = convnfft_fast(squeeze(image), gaussian_kernel(0.65));
warning('CellOrganizer:getMeshPoints', 'Smoothing temporarily disabled due to memory issues, see utilities/3D/vesicles/3D/tp_imresize3d.m for example of estimating memory usage');
smoothing_time = toc(smoothing_time);

% Downsampling
disp('Downsampling image')
downsampling_time = tic;
% Runs out of memory on lanec1 model1
% D = tp_imresize3d(image,round(size(image)*downsampling),'linear');
warning('CellOrganizer:getMeshPoints', 'Downsampling temporarily using nearest due to memory issues');
D = tp_imresize3d(image,round(size(image).*downsampling),'nearest');
downsampling_time = toc(downsampling_time);

%make iso-surface (Mesh) of skin
save([options.temporary_results,filesep,'D.mat'],'D');
disp('Creating mesh')
isosurface_time = tic;
FV = isosurface(D,isopatch);
isosurface_time = toc(isosurface_time);
disp('Adjusting mesh')
%{
%Shift the vertices back to where they're supposed to be so they line
%up with the objects inside the cell.
%add the amount of shift, and subtract the padding
%get the top value of the real image
%Must subtract 1 from them to create zero indexing
topval = [max(xvals)-1,max(yvals)-1,max(zvals)-1];
shift = round(max(FV.vertices)-topval);
%}
% Shift to original coordinates, assumes all images are the same size
shift = padAmount-[min(xvals),min(yvals),min(zvals)];
% This is close but not exact because image size is rounded for tp_imresize3d
shift = shift.*downsampling;
FV.vertices(:,1) = FV.vertices(:,1)-shift(1);
FV.vertices(:,2) = FV.vertices(:,2)-shift(2);
FV.vertices(:,3) = FV.vertices(:,3)-shift(3);


% Compensate for downsampling
% FV.vertices = FV.vertices ./ downsampling;
for i = 1:3
    FV.vertices(:,i) = FV.vertices(:,i) ./ downsampling(i);
end


% FV.vertices = FV.vertices-repmat(shiftvector_flag,size(FV.vertices,1),1);
fprintf('Before reducepatch: %i vertices, %i faces\n', size(FV.vertices, 1), size(FV.faces, 1));
FV = reducepatch(FV,patchsample);
fprintf('After reducepatch: %i vertices, %i faces\n', size(FV.vertices, 1), size(FV.faces, 1));
