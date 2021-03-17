function answer = demo3D59( options )
% demo3D59
%
% Synthesize one 3D image with nuclear, cell shape and a vesicular channel.
% This demo exports portions of the synthetic image as an SBML Spatial instance.
%
% Input
% -----
% * a valid CellOrganizer model
%
% Output
% ------
% * SBML instance
% * single channel TIF files

% Author: Ivan E. Cao-Berg, Taraz Buck
%
% Copyright (C) 2016-2019 Murphy Lab
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
% FEEL FREE TO MODIFY THE VARIABLES IN THIS BLOCK

% The framework (cell and nucleus shapes) will be synthesized using framework_model
framework_model = '../../../models/3D/spharm/lamp2.mat';
vesicle_models = '../../../models/3D/tfr.mat';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT MODIFY THIS BLOCK
if ~isdeployed()
  current_path = which(mfilename);
  [current_path, filename, extension] = fileparts( current_path );
  cd(current_path);
end

disp( mfilename );
disp( 'The estimated running time is about 6 minutes. Please wait.' );

% Combine SPHARM framework with vesicle models
combined_models = 'combined_model.mat';
slml2slml({framework_model, vesicle_models}, struct('output_filename', combined_models, 'selection', [1,1,0;0,0,1]));

options.seed = 639848;
options.targetDirectory = pwd;
options.prefix = 'img';
options.synthesis = 'all';
options.model.spharm_rpdm.synthesis_method = 'random_sampling';
options.model.spharm_rpdm.imageSize = [205, 205, 18];
options.numberOfSynthesizedImages = 1;
options.numberOfGaussianObjects = 100;
options.rendAtStd = 1;
options.overlapthresh = 1;
options.compression = 'lzw';
options.microscope = 'none';
options.sampling.method = 'disc';
options.verbose = true;
options.debug = false;
options.output.tifimages = true;
options.output.OMETIFF = true;
options.output.indexedimage = true;
options.output.SBMLSpatial = true;
% options.SBML.downsampling = [1/4, 1/4, 1];
options.SBML.spatial_use_compression = true;
options.SBML.spatial_use_analytic_meshes = ~false;
options.SBML.spatial_image = false;
options.SBML.spatial_vcell_compatible = false;
options.oobbuffer = 0.1;

try
    state = rng( options.seed );
catch err
    rand( 'seed', options.seed ); %#ok<RAND>
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

answer = slml2img( {combined_models}, options );
end%demo3D59
