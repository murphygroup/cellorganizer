function answer = demo3D17()
% demo3D17
%
% The main idea behind this demo is to show the user they
% can use their own binary images from raw experimental data 
% to synthesize protein patterns. 
% 
% The current demo assumes the resolution of the images is the same 
% as the resolution of the images that were used to train the protein model.
%
% Input 
% -----
% * an existing raw or synthetic framework, i.e. one binary multi-TIFF
% file of the nuclear channel and one binary multi-TIFF file of the
% cell membrane
% * the resolution of the latter images
% * a valid CellOrganizer model that contains a protein model
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

disp( 'demo3D17' );
disp( 'The estimated running time is 2 minutes. Please wait.' );

try
    state = rng(3);
catch
    state = RandStream.create('mt19937ar','seed',3);
    RandStream.setDefaultStream(state);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FEEL FREE TO MODIFY THE VARIABLES IN THIS BLOCK
%step0: set the input parameters described in the header
%step0.1 use existing raw or synthetic frameworks: in this demo we are going
%to use an image from the Murphy Lab 3D HeLa collection that you can download
%along with CellOrganizer.

filename = '../../../images/HeLa/3D/processed/LAM_cell10_ch0_t1.tif';
param.instance.nucleus = tif2img( filename );
filename = '../../../images/HeLa/3D/processed/LAM_cell10_ch1_t1.tif';
param.instance.cell = tif2img( filename );

%step0.2: set the resolution of the latter images
param.instance.resolution = [0.049, 0.049, 0.2000];

%step0.3: use a valid CellOrganizer model that contains a protein model. in
%this model we are going to use the 3D HeLa nucleoli model distrubuted in
%this version of CellOrganizer
model_file_path = '../../../models/3D/nuc.mat';

%these are optional parameters that you are welcome to modify as needed
%location where CellOrganizer will save the images to
param.targetDirectory = pwd;
param.prefix = 'img';
param.numberOfSynthesizedImages = 1;
param.output.tifimages = true;
param.compression = 'lzw';
param.microscope = 'none';
param.sampling.method = 'disc';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

answer = slml2img( {model_file_path},  param );
