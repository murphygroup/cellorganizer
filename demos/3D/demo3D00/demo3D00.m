function answer = demo3D00()
% demo3D00
%
% Synthesize one 3D image with nuclear, cell shape, and nucleolar channels
% from nucleolar model with sampling method set to render nucleoli as
% ellipsoids without convolution. The model was trained from the Murphy Lab
% 3D HeLa dataset.
%
% Input
% -----
% * a valid CellOrganizer model file
%
% Output
% ------
% * three TIFF files (nuclear, cell shape, and nucleolar channels)

% Robert F. Murphy
%
% Copyright (C) 2012-2019 Murphy Lab
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

disp( 'demo3D00' );
disp( 'The estimated running time is 1 minute. Please wait.' );

options.seed = 12345;
try
    state = rng( options.seed );
catch err
    rand( 'seed', options.seed ); %#ok<RAND>
end

options.targetDirectory = pwd;
options.prefix = 'img';
options.compression = 'lzw';
options.sampling.method = 'disc';
options.microscope = 'none';
options.debug = false;
options.verbose = true;
options.display = false;
options.output.OMETIFF = true;
options.output.tifimages = true;
options.temporary_results = [ pwd filesep 'temporary_results' ];
options.overwrite_synthetic_instances = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FEEL FREE TO MODIFY THE VARIABLES IN THIS BLOCK
list_of_models = {'../../../models/3D/nuc.mat'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

answer = slml2img( list_of_models, options );
end%demo3D00
