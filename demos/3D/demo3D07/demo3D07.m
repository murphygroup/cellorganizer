function answer = demo3D07()
% demo3D07
%
% Synthesize one 3D image with nuclear, cell shape, and protein channels
% from all object models (nucleoli, lysosomes, endosomes, mitochondria, and
% microtubules) with sampling method set to sample vesicular objects from
% Gaussians at a density of 25 and convolution with point-spread function.
% The model was trained from the Murphy Lab 3D HeLa dataset.
%
% Input 
% -----
% * a list of valid CellOrganizer model files
%
% Output
% ------
% * seven TIFF files (nuclear, cell shape, nucleolar, lysosomal, endosomal,
%   mitochondrial, and microtubule channels)

% Robert F. Murphy
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

disp( 'demo3D07' );
disp( 'The estimated running time is 30 sec. Please wait.' );

options.seed = 100;
try
    state = rng( options.seed );
catch err
    rand( 'seed', options.seed ); %#ok<RAND>
end

options.numberOfSynthesizedImages = 1;
options.targetDirectory = pwd;
options.prefix = 'img';
options.compression = 'lzw';
options.microscope = 'svi';
options.sampling.method = 'sampled';
options.sampling.density = 25;
options.verbose = true;
options.debug = false;
options.numberOfGaussianObjects = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FEEL FREE TO MODIFY THE VARIABLES IN THIS BLOCK
list_of_models = {'../../../models/3D/nuc.mat', ...
    '../../../models/3D/lamp2.mat', ...
    '../../../models/3D/tfr.mat', ...
    '../../../models/3D/mit.mat', ...
    '../../../models/3D/centro.mat'...
    '../../../models/3D/tub.mat'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

answer = slml2img( list_of_models, options );
