function [indexedimg,compartmentlist] = DefineCompartments(cellimg,nucimg,param)
%This function uses a cell and nuclear segmentation to define cellular
%compartment masks. the result is an indexed image where each index
%corresponds to a compartment. 
%
%Inputs: 
%
%cellimg = image of segmented cell (2 or 3D)
%nucimg = image of segmented nuclei(2 or 3D)
%param = struct array containing the resolution of the images 
%
%Outputs: 
%
%indexedimg = resulting image where each integer pixel value corresponds to
%             a different compartment
%compartmentlist = a string array of compartments in order

%Author: Devin Sullivan 5/30/13
% G. Johnson 7/7/13 - Improved membrane detection speed by 100x by using
% bwdist
%
% Copyright (C) 2012 Murphy Lab
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

%initialize image
indexedimg = zeros(size(cellimg));

%this is hard coded now, but the idea is that in the future, the user could
%specify which compartments he/she wants and it will only grab those
currentidx = 1;

%%%CYTO%%%
indexedimg(cellimg>0)=currentidx;
compartmentlist{currentidx} = 'cyto';
currentidx = currentidx+1;


%%%NUCLEUS%%%
indexedimg(nucimg>0) = currentidx;
compartmentlist{currentidx} = 'nuc';
currentidx = currentidx+1;


%%%MEMBRANE%%%
%This is the trickiest one as the width you wish to consider depends on
%your resolution. Here we will make the assumption that the plasma membrane
%is any area within 0.1um of the cell edge. 

% %First get the cell edge 
% celledge = bwperim(cellimg>0);
% %then get the xyz positions of each point 
% [xpos,ypos,zpos] = ind2sub(size(cellimg),[1:length(cellimg(:))]);
% %for some reason xyz needs to be transposed because it comes out sideways
% xyz = [xpos;ypos;zpos]';
% 
% %get positions of only values that are on the perimeter
% [xp,yp,zp] = ind2sub(size(celledge>0),find(celledge(:)>0));
% xyz_perim = [xp,yp,zp];


%now find all the points within a given distance (we set this using the
%resolution, but in the future the user could set it)
%set distance to the floor in pixels equal to 0.01um. This value is chosen
%from literature as the approximate thickness of the plasmamembrane.
% memthick = 0.01;%um
%however we would like to capture the formation of vesicles which could be
%substantially larger than the membrane thickness, so we introduce a
%parameter for the size of the vesicles. In the future, this could be taken
%from the cell specific model. Here we take the number from Stoorvogel et
%al JCB 1996.  of 100nm=0.1um
pitthick = 0.1;%um
size_membrane = pitthick./param.model.resolution;

%(0.1um)/(Xum/pixel)= pixels/0.1um, allow for different resolutions in each
%dimension

% %find indices within the distance (note, this is the slow step)
% ind = rangesearch(xyz,xyz_perim,1,'Distance','seuclidean','Scale',size_membrane);
% 
% %loop through the ind cellarray and concatinate them 
% totind = [];
% for i = 1:length(ind)
%     totind = [totind,ind{i}];
% end

%assign values to make the mask
%don't actually need the membrane image, but could be useful in the future
% cellmem = zeros(size(cellimg));
% cellmem(totind) = 1;
% indexedimg(totind) = currentidx;
% compartmentlist{currentidx} = 'mem';

%grj 7/7/13 - Improved speed over rangesearch by 100x
dist = bwdist(~padarray(cellimg,[0,0,1], 0, 'both'));
dist = dist(:,:,2:end-1);

indexedimg(dist <= ceil(size_membrane(1)) & cellimg) = currentidx;
compartmentlist{currentidx} = 'mem';
currentidx = currentidx+1;






%in the future it would be interesting to see how much of the current
%pattern fluorescence is fit by the gaussian objects 
%here I outline code to do this for future use
%**NOTE: This will overwrite current mask indexed images (e.g. cytoplasmic
%space).If just testing amount of fluor in pattern use a separate image
%this case is dealt with in CompartmentProportions directly, however it
%does break the dirchilet (summing to 1)
%if param.protmask==1
% indexedimg(param.protmask>0) = currentidx;
% compartmentlist{currentidx} = 'prot';
% currentidx = currentidx+1;
%end


