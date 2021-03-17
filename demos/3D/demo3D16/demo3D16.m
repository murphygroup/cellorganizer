function answer = demo3D16()
% demo3D16
%
% The main idea behind this demo is to show the user they
% can use their own binary images from raw experimental data 
% to synthesize protein patterns. This demo uses the CellOrganizer
%  method for nuclear and cell segmentation.
% 
% The current demo assumes the resolution of the images is the same as 
% the resolution of the images that were used to train the protein model.
%
% Input 
% -----
% * raw or synthetic images of the nuclear and cell membrane
% * a valid CellOrganizer model file
%
% Output
% ------
% * three TIFF files (cell shape, nuclear, and lysosomal channels)

% Ivan E. Cao-Berg
%
% Copyright (C) 2012-2017 Murphy Lab
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
  cd(current_path);
end

disp( 'demo3D16' );
disp( 'The estimated running time is 3 hours. Please wait.' );

options.seed = 12345;
try
    state = rng( options.seed );
catch err
    rand( 'seed', options.seed ); %#ok<RAND>
end

options.targetDirectory = pwd;
options.prefix = 'img';
options.numberOfSynthesizedImages = 1;
options.compression = 'lzw';
options.microscope = 'none';
options.sampling.method = 'disc';
options.verbose = true;
options.debug = false;
options.display = false;
options.synthesis = 'all';
options.train.flag = 'all';
options.resolution.cell = [0.049, 0.049, 0.2000];
downsample = [1 1 1];
options.resolution.cell = options.resolution.cell.*downsample;
options.preprocessingFolder = 'temp';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

directory = '../../../images/HeLa/3D/raw';
dna_image_file = [ directory filesep 'LAM_cell1_ch0_t1.tif' ];
cell_image_file = [ directory filesep 'LAM_cell1_ch1_t1.tif' ];
prot_image_file = [ directory filesep 'LAM_cell1_ch2_t1.tif' ];
mask_image_file = [ directory filesep 'LAM_cell1_mask_t1.tif' ];
options.instance.resolution=[0.049, 0.049, 0.2000];

[options.instance.nucleus, options.instance.cell ] = ...
    seg_cell_and_nucleus( dna_image_file, ...
    cell_image_file, ...
    prot_image_file, ...
    mask_image_file, ...
    downsample, options.display, options );

answer = slml2img( {'../../../models/3D/lamp2.mat'}, options );
