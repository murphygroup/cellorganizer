function answer = demo3D15()
% demo3D15
%
% Synthesize one multichannel 3D image from an endosomal model and
% diffeomorphic nuclear and cell shape model. The sampling method was set
% to render endosomes as ellipsoids without convolution. The model was
% trained from the Murphy Lab 3D HeLa dataset.
%
% Input 
% -----
% * a valid CellOrganizer model file with a diffeomorphic framework
%
% Output
% ------
% * three TIFF files (nuclear, cell shape, and endosomal channels)

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

disp( 'demo3D15' );
warning('This demo is deprecated. The demo will be removed in future versions of CellOrganizer');

options.seed = 1234558;
try
    state = rng( options.seed );
catch err
    rand( 'seed', options.seed ); %#ok<RAND>
end

outputDirectory = pwd;
options = [];
options.targetDirectory = outputDirectory;
options.prefix = 'img';
options.numberOfSynthesizedImages = 1;
options.numberOfGaussianObjects = 5;
options.compression = 'lzw';
options.microscope = 'none';
options.sampling.method = 'disc';
options.verbose = true;
options.debug = false;
options.synthesis = 'all';

filename = '../../../models/3D/diffeomorphic/HeLa_diffeo_framework_and_LAMP2.mat';
if ~exist( filename )
    warning(['Unable to load model: ' filename ' . Exiting method.'] );
    return
else
    load( filename );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FEEL FREE TO MODIFY THE VARIABLES IN THIS BLOCK
list_of_models = { filename };
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

answer = slml2img( list_of_models, options );
