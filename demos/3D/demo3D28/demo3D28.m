function answer = demo3D28()
% demo3D28
%
% Synthesize one 3D image with nuclear, cell shape, and nucleolar channels
% from nucleolar model with sampling method set to render nucleoli as
% ellipsoids without convolution. The model was trained from the Murphy Lab
% 3D HeLa dataset.
%
% Input
% -----
% * an existing raw or synthetic nuclear image, i.e. one binary multi-TIFF
%   file of the nuclear channel
% * the resolution of the input image
% * a valid CellOrganizer model that contains a cell membrane model
%
% Output
% ------
% * three TIFF files (cell shape, nuclear, and nucleolar channels)

% Ivan E. Cao-Berg
%
% Copyright (C) 2016-2017 Murphy Lab
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

disp( 'demo3D28' );
disp( 'The estimated running time is 1 minute. Please wait.' );

options.seed = 3;
try
    state = rng(options.seed);
catch
    state = RandStream.create('mt19937ar','seed',options.seed);
    RandStream.setGlobalStream(state);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FEEL FREE TO MODIFY THE VARIABLES IN THIS BLOCK
%step0: set the input parameters described in the header
%step0.1 use existing raw or synthetic frameworks: in this demo we are going
%to use an image from the Murphy Lab 3D HeLa collection that you can download
%along with CellOrganizer.
filename = '../../../images/HeLa/3D/processed/LAM_cell10_ch0_t1.tif';
options.instance.nucleus = tif2img( filename );

%step0.2: set the resolution of the latter images
options.instance.resolution = [0.049, 0.049, 0.2000];

%step0.3: use a valid CellOrganizer model that contains a protein model. in
%this model we are going to use the 3D HeLa nucleoli model distributed in
%this version of CellOrganizer
model_file_path = '../../../models/3D/nuc.mat';

%these are optional parameters that you are welcome to modify as needed
%location where CellOrganizer will save the images to
options.targetDirectory = pwd;

%output folder name
options.prefix = 'img';

%number of images to synthesize
options.numberOfSynthesizedImages = 1;

%save images as TIF files
options.output.tifimages = true;

%compression for TIF output
options.compression = 'lzw';

%do not apply point-spread-function
options.microscope = 'none';

%render Gaussian objects as discs
options.sampling.method = 'disc';

%synthesize for all channels
options.synthesis = 'all';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%main call to CellOrganizer
answer = slml2img( {model_file_path},  options );
end%demo3D28