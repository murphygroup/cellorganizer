function answer = demo3D31()
% demo3D31
%
% Trains a generative model of microtubules
%
% Input 
% -----
% * a directory of raw or synthetic nucleus images
% * a directory of raw or synthetic cell shape images
% * the resolution of the images (all images should have the same
%   resolution)
%
% Output
% ------
% * a valid model

% Ivan E. Cao-Berg & Gregory R. Johnson
%
% Copyright (C) 2014-2017 Murphy Lab
% Computational Biology Department
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT MODIFY THIS BLOCK
if ~isdeployed()
  current_path = which(mfilename);
  [current_path, filename, extension] = fileparts( current_path );
  get_murphylab_image_collections( true );
  cd(current_path);
end

disp( 'demo3D31' );
disp( 'The estimated running time is greater than 12 hours. Please be patient.' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pattern = 'Tub';
imgdir = '../../../images/HeLa/3D/raw';
options.protein.class = 'network';
options.protein.type = 'microtubule_growth';
options.nucleus.class = 'nuclear_membrane';
options.nucleus.type = 'cylindrical_surface';
options.cell.class = 'cell_membrane';
options.cell.type = 'ratio';
options.protein.name = 'Tub';
options.nucleus.name = 'Tub';
options.cell.model = 'Tub';
options.faster = false;
options.train.flag = 'all';

dimensionality = '3D';
dna = [imgdir filesep 'Tub*cell*ch0*.tif'] ;
cell = [imgdir filesep 'Tub*cell*ch1*.tif'];
protein = [imgdir filesep 'Tub*cell*ch2*.tif'];
mask = [imgdir filesep 'Tub*cell*mask*.tif'];

% generic model options
% ---------------------
options.model.name = 'microtubule_growth';
options.model.id = num2str(now);
options.model.filename = [ lower(pattern) '.xml' ];

% nuclear shape model options
% ---------------------------
options.nucleus.class = 'nuclear_membrane';
options.nucleus.type = 'cylindrical_surface';
options.nucleus.id = num2str(now);

% cell shape model options
% ------------------------
options.cell.class = 'cell_membrane';
options.cell.type = 'ratio';
options.cell.id = num2str(now);

% protein shape model options
% ---------------------------
options.protein.class = 'network';
options.protein.type = 'microtubule_growth';
options.protein.id = num2str(now);

% other options
% -------------
options.verbose = true;

options.debug = true;
%options.train.flag = 'microtubule';
options.model.resolution = [0.049, 0.049, 0.2000];
options.downsampling = [4,4,1];
options.saveIntermediateResults = false;
options.model.microtubule.searchparams.n = [50:100:500];
options.model.microtubule.searchparams.mulen = [10:10:50];
options.model.microtubule.searchparams.colli_min_number = [0.97, 1 ];
options.protein.cytonuclearflag = 'cyto';
options.debug = true;
options.verbose = true;

dna_paths = dir(dna);
dna_paths = sort_nat({dna_paths.name});

prot_paths = dir(protein);
prot_paths = sort_nat({prot_paths.name});

cell_paths = dir(cell);
cell_paths = sort_nat({cell_paths.name});

mask_paths = dir(mask);
mask_paths = sort_nat({mask_paths.name});

for i = 1:length(dna_paths)
    cell_imgs{i} = @() seg_hela([imgdir filesep cell_paths{i}], [imgdir filesep mask_paths{i}]);
    dna_imgs{i} = @() seg_hela([imgdir filesep dna_paths{i}], [imgdir filesep mask_paths{i}]);
    prot_imgs{i} = [imgdir filesep prot_paths{i}];
end

if ~exist([pwd filesep 'tub.mat'], 'file')
    answer = img2slml( dimensionality, dna_imgs, cell_imgs, prot_imgs, options );
else
    answer = true;
end
end%demo3D31

function [ seg_region ] = seg_hela( impath, maskpath )
%SEG_HELA Summary of this function goes here
%   Detailed explanation goes here
img = ml_readimage(impath);
mask = ml_readimage(maskpath);
mask = mask(:,:,1);
mask = repmat(mask, [1,1, size(img,3)]) > 0;
img = mat2gray(img);
img_region = img > graythresh(img);
img_region = img_region.*mask;
regions = bwconncomp(img_region);
[~, ind] = max(cellfun(@(x) size(x,1), regions.PixelIdxList));
seg_region = false(size(img));
seg_region(regions.PixelIdxList{ind}) = true;
seg_region = imfill(seg_region, 'holes');
end%seg_hela
