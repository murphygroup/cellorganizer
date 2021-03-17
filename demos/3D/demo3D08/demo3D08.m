function answer = demo3D08()
% demo3D08
%
% Synthesize one 3D image with nuclear, cell shape, and vesicular channels
% from all vesicular object models (nucleoli, lysosomes, endosomes, and
% mitochondria) with sampling method set to render vesicular objects as
% ellipsoids without convolution. The model was trained from the Murphy Lab
% 3D HeLa dataset.
%
% Input 
% -----
% * a list of valid CellOrganizer model files
%
% Output
% ------
% * single indexed TIFF file which indexes the six TIFF files (nuclear,
%   cell shape, nucleolar, lysosomal, endosomal, and mitochondrial channels)

% Created: Michelle Mackie
%
% Copyright (C) 2012-2017 Murphy Lab
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT MODIFY THIS BLOCK
if ~isdeployed()
  current_path = which(mfilename);
  [current_path, filename, extension] = fileparts( current_path );
  cd(current_path);
end

disp( 'demo3D08' );
disp( 'The estimated running time is 4 minutes. Please wait.' );

options.seed = 100;
try
    state = rng( options.seed );
catch err
    rand( 'seed', options.seed ); %#ok<RAND>
end

options.numberOfSynthesizedImages = 1;
options.targetDirectory = pwd;
options.prefix = 'img';
options.image.compression = 'lzw';
options.microscope = 'none';
options.sampling.method = 'disc';
options.verbose = true;
options.debug = false;
options.output.indexedimage = true;
options.overlapthresh = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FEEL FREE TO MODIFY THE VARIABLES IN THIS BLOCK
list_of_models = {'../../../models/3D/nuc.mat', ...
    '../../../models/3D/lamp2.mat', ...
    '../../../models/3D/tfr.mat', ...
    '../../../models/3D/mit.mat'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

answer = slml2img( list_of_models, options );
