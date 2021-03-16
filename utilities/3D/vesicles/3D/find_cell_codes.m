function [cellcodes, equatorZ, cellnucheightratio, nucbottomslice] = find_cell_codes( seg_dna, seg_cell, param )
%FIND_CELL_CODES Calculates cell codes

% Tao Peng
%
% Copyright (C) 2011-2016 Murphy Lab
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

% March 12, 2012 Added another preprocessing routine. This routine mimics
%                the preprocessing from A. Shariff, G. K. Rohde and R. F. Murphy (2010) A
%                Generative Model of Microtubule Distributions, and Indirect Estimation of
%                its Parameters from Fluorescence Microscopy Images. Cytometry Part A 77A:457-466.
%
% March 19, 2012 I. Cao-Berg Added Active Countour 3D segmentation
%
% March 22, 2012 I. Cao-Berg Added param.downsample so that user can select a
%                downsampling scale. The default is [1 1 1] which is
%                the one in the original generative model paper
%
% March 23, 2012 I. Cao-Berg Fixed a bug when getting cellnum from impath
%
% March 27, 2012 I. Cao-Berg When calling ml_parsecell instead of sending a slice
%                of the nucleus, now we are sending a z-projection of the nucleus
%
% April 5, 2012 I. Cao-Berg Added new preprocessing method
%
% April 9, 2012 R.F. Murphy Save cell-nucleus height ratio and nuc bottomslice
%
% June 19, 2012 I. Cao-Berg Changed index name in the preprocessing for-loops that had the same
%                           index name as the main loop
%
% July 26, 2012 I. Cao-Berg If the number of slices in the nuclear image that have
%                           fluorescence is less than 4, then ignore such image
%
% May 15, 2013 I. Cao-Berg Updated method to support wildcards
%
% June 4, 2013 D. Sullivan Added check if masks were specified to avoid
%                          index out of bounds code when they're not
%
%
% June 7-13 2013 D. Sullivan changed to per-cell type computations. Major
%                            changes include passing single files instead
%                            of directories and passing in a tempdirectory
%                            structure as part of param.
%
% March 3, 2016 I. Cao-Berg Updated default flags' values

if nargin == 2
    param = [];
end

param = ml_initparam( param, struct( 'display', false ) );
param = ml_initparam( param, struct( 'debug', false ) );
param = ml_initparam( param, struct( 'verbose', false ) );

%create empty containers
cellcodes = [];
equatorZ = [];
cellnucheightratio = []; 
nucbottomslice = [];

nucbody = {};
cellbody = {};

%icaoberg 26/7/2012
nnzpixel = squeeze(sum(sum(seg_dna,1),2));

if length(find(nnzpixel)>0) < 4
    if param.verbose
        disp(['Ignoring image']);
    end
    %D. Sullivan 6/6/13 want to return instead of continue now
    %that this code does one cell at a time.
    %                 continue
    return
end


%icaoberg april 5, 2012
stacknumber = 15;
[seg_dna, seg_cell] = ml_rescaleImage2Cell( seg_dna(:,:,nnzpixel>0), seg_cell(:,:,nnzpixel>0), stacknumber );

%fills the empty slices of DNA with the closest
z = find(sum(sum(seg_dna,1),2));
seg_dna(:,:,1:min(z)-1)=repmat(seg_dna(:,:,min(z)),[1,1,min(z)-1]);
seg_dna(:,:,max(z)+1:end)=repmat(seg_dna(:,:,max(z)),[1,1,size(seg_dna,3)-max(z)]);

%needed for model of nuclear position
cellnucheightratio = stacknumber/length(z);
nucbottomslice = (min(z)-1)/stacknumber;

%preprocess DNA image
se = strel('disk',3);
crossectionArea = [];
%icaoberg june 19, 2012 changed index from i to j
for j=1:1:size(seg_dna,3)
    mask = seg_dna(:,:,j);
    mask = imclose(mask,se);
    mask = imfill(mask,'holes');
    obj = ml_findmainobj2d_bw(mask);
    crossectionArea(j) = nnz(mask);
    nucbody{j}=obj;
end
[maxArea,equatorZ] = max(crossectionArea);

%preprocess cell image
%icaoberg june 19, 2012 changed index from i to j
for j=1:1:size(seg_cell,3)
    mask = seg_cell(:,:,j);
    mask = imclose(mask,se);
    mask = imfill(mask,'holes');
    if ~all(mask(:) == 0)
        obj = ml_findmainobj2d_bw(mask);
        cellbody{j}=obj;
    end
    
end

cellcodes = {};
for k = 1:1:length(cellbody)
    if param.verbose
        disp(['Calculating cell code on stack number: ' num2str(k) ]);
    end
    
    method = 0;
    if method == 0
        %calculate projection here
        %let nucbody become a z-projection and calculate distance according to the projection
        mask = sum( seg_dna, 3 );
        mask(find( mask ~= 0)) = 255;
        mask = logical( mask );
        obj = ml_findmainobj2d_bw(mask);
        crossectionArea = nnz(mask);
        nucbody = obj;
        
        if ~isempty(nucbody) && ~isempty(cellbody{k})
            cellcodes{k} = ml_parsecell({},cellbody{k},nucbody,...
                1,size(seg_dna(:,:,k)),...
                {'da','nucarea','nuccenter','nucmangle','nuchitpts',...
                'nuccontour','nucellhitpts','nucdist','nucelldist','nucecc',...
                'cellarea','cellcenter','cellmangle','cellcontour',...
                'cellhitpts','celldist','cellecc','nucedge','celledge'},0);
        end
    elseif method == 1
        %calculate projection here
        %let nucbody become a z-projection and calculate distance according to the projection
        %                     mask = sum( seg_dna, 3 );
        %                     mask(find( mask ~= 0)) = 255;
        %                     mask = logical( mask );
        %                     obj = ml_findmainobj2d_bw(mask);
        %                     crossectionArea = nnz(mask);
        %nucbody = obj;
        nuclearbody = nucbody{k};
        
        if ~isempty(nucbody)
            cellcodes{k} = ml_parsecell({},cellbody{k},nuclearbody,...
                1,size(seg_dna(:,:,k)),...
                {'da','nucarea','nuccenter','nucmangle','nuchitpts',...
                'nuccontour','nucellhitpts','nucdist','nucelldist','nucecc',...
                'cellarea','cellcenter','cellmangle','cellcontour',...
                'cellhitpts','celldist','cellecc','nucedge','celledge'},0);
        end
    elseif method == 2
        %calculate projection here
        %let nucbody become a z-projection and calculate distance according to the projection
        mask = sum( seg_dna, 3 );
        mask(find( mask ~= 0)) = 255;
        mask = logical( mask );
        obj = ml_findmainobj2d_bw(mask);
        crossectionArea = nnz(mask);
        %nucbody = obj;
        nuclearbody = nucbody{k};
        
        if ~isempty(nucbody)
            cellcodes{k} = ml_parsecell({},cellbody{k},nuclearbody,...
                1,size(seg_dna(:,:,k)),...
                {'da','nucarea','nuccenter','nucmangle','nuchitpts',...
                'nuccontour','nucellhitpts','nucdist','nucelldist','nucecc',...
                'cellarea','cellcenter','cellmangle','cellcontour',...
                'cellhitpts','celldist','cellecc','nucedge','celledge'},0);
        end
    end
end
end

