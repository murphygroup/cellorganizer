function [param] = findobjs_run(...
    dnaImagesDirectoryPath, ...
    cellImagesDirectoryPath, ...
    protImagesDirectoryPath,savepath,param)
% Extract objects from 3D HeLa images

% Tao Peng
%
% Copyright (C) 2011-2013 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
% Carnegie Mellon University
%
% June 15, 2012 G. Johnson Changed ml_findobjs loop to accept images of
%         arbritrary x and y dimension
% July 18, 2012 I. Cao-Berg Updated method to ignore single pixel objects
% August 1, 2012 I. Cao-Berg Fixed a bug in the code that would insert an artifact to the 
%                            end of the image that would later be considered by the 
%                            algorithm as an object
% August 2, 2012 D. Sullivan Added masking of protein image using segcell.
%                            this ensures you find no objects outside of
%                            the segmented cell
% Feb 22, 2013 D. Sullivan   Added protein resolution adjustment so that
%                            gaussians are trained on cubic voxels with the
%                            maximum resolution
% Feb 24, 2013 D. Sullivan   Removed Feb 22 change for speed
% May 15, 2013 I. Cao-Berg Updated method to support wildcards
%
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

if ~exist(savepath,'dir')
    mkdir(savepath)
end

%icaoberg 15/5/2013
output_directory = [savepath filesep 'preprocessed'];
if ~exist( output_directory, 'dir' )
    mkdir( output_directory );
end

%icaoberg 15/5/2013
dna_image_files = ml_ls( dnaImagesDirectoryPath );
cell_image_files = ml_ls( cellImagesDirectoryPath );
prot_image_files = ml_ls( protImagesDirectoryPath );
try
    masks_image_files = ml_ls( param.masks );
catch
    mask_image_files = '';
end

H = fspecial3('Gaussian',[10,10,3]);

for i = 1:length(prot_image_files)
    disp(['Image:' num2str(i)])
    %icaoberg 15/5/2013
    try
        dna_image_file = dna_image_files{i};
    catch
        dna_image_file = '';
    end
    
    cell_image_file = cell_image_files{i};
    prot_image_file = prot_image_files{i};
    
    try
%         mask_image_file = masks_image_files{i};
        crop = ml_loadimage(mask_image_files{i});
    catch
%         mask_image_file = [];
        crop = [];
    end
    
    if ~exist( [savepath filesep 'obj' num2str(i) '.mat'], 'file' )
     try
        prot = ml_readimage( prot_image_file );
        stacknum = size(prot,3);
        prot = double(prot);
        
        %%%D. Sullivan 2/22/13 
        %Changing how training works. Resolution is now used to specify the
        %model
%         %make sure the protein pattern has the same resolution in z as it does in x&y
%         xsize = floor(param.model.original_resolution(1)/min(...
%             param.model.original_resolution)*size(prot,1));
%         ysize = floor(param.model.original_resolution(2)/min(...
%             param.model.original_resolution)*size(prot,2));
%         zsize = floor(param.model.original_resolution(3)/min(...
%             param.model.original_resolution)*size(prot,3));
%         prot = imresize(prot,[xsize,ysize]);
%         prot = tp_stretch3d(prot,zsize);
%         param.model.protein_resolution(1:3) = min(param.model.original_resolution);
        %%%
        
        %D. Sullivan 6/4/13 get rid of this due to new file formats (bioformats) 
%         croplist = ml_ls( mask_image_file );
%         crop = [];
%         for k = 1:length(croplist)
%             crop(:,:,k) = ml_readimage(croplist{k}) > 0;
%         end
        crop = any(crop,3);
        
        disp( 'Preprocessing image' )
        %devins 8/2/12
        %Need to make sure we are not finding objects outside of the cell
        segdataFolder = [ pwd filesep 'temp' filesep 'preprocessing' filesep 'cell' num2str(i) '.mat' ];

        %this file should contain variables 'segcell', 'segdna', 'downsample'
        load(segdataFolder);
        %need to make the image it's original size
        %create a blank image
        cellmask = zeros(original_segcellsize);
        %put segcell into the blank image

        cellmask(:,:,bot_slice:top_slice) = segcell;
        %figure out what size to make the final image
        finalsize = size(prot);
        %resize x and y
        resizedcell= imresize(cellmask,finalsize(1:2));
        %shouldn't need to resize the z since it does not change in
        %preprocessing, but just in case...
        %make sure the z dimension is the correct size 
        resizedcell = tp_stretch3d(resizedcell,finalsize(3));
        
        %create the masked protein image
        prot = prot.*logical(resizedcell);
        
        %end of 8/2/12 addition
        %%%
        
        procimage = ml_preprocess(double(prot),crop,'ml','yesbgsub','nih');
        %procimage(1:12,:,:) = [];
        %procimage(end-11:end,:,:) = [];
        %procimage(:,1:12,:) = [];
        %procimage(:,end-11:end,:) = [];
        
        dimobjimg = imfilter(procimage,H,'same');
        centers{i} =  ml_imlocalmax(dimobjimg);
        
        save([output_directory filesep 'cell' num2str(i) '.mat'], ...
               'procimage', 'prot', 'crop', 'resizedcell', ...
               'dna_image_file', 'cell_image_file', 'prot_image_file', ...
               'mask_image_file' );

        blockSize = 100;
        numBlocksX = round(size(procimage,1) / blockSize);
        numBlocksY = round(size(procimage,2) / blockSize);
        
        objects = cell(numBlocksX, numBlocksY);
        for r = 1:numBlocksX
            for c = 1:numBlocksY
                disp( ['Finding objects in block {' num2str(r) ',' num2str(c) '}' ] );
                if r == numBlocksX
                    %icaoberg 1/8/2012
                    rEnd = size(procimage,1)-1;
                else
                    rEnd = r*blockSize;
                end
                if c == numBlocksY
                    %icaoberg 1/8/2012
                    cEnd = size(procimage,2)-1;
                else
                    cEnd = c*blockSize;
                end
              
                %icaoberg 17/7/2012
                blockimage = procimage((r-1)*blockSize+1:rEnd,(c-1)*blockSize+1:cEnd,:);
                objs = ml_findobjs(blockimage);
                objSizes = cellfun('size',objs,1);
                objs = objs(objSizes>1);
                objects{r,c} = objs;
            end
        end
        
        save([savepath filesep 'obj' num2str(i)],'objects','centers', 'blockSize', 'numBlocksX', 'numBlocksY')
     catch
        disp(['Unable to process image ' prot_image_file ]);
     end
    else
        disp( 'Intermediate files found, skipping locating objects' );
    end
end
