function answer = demo3D60( options )
% demo3D60
%
% Synthesize one 3D image with nuclear, cell shape and a vesicular channel using the SPHARM framework model.
% This demo exports the synthetic image as Virtual Cell VCML.
%
% Input 
% -----
% * a valid CellOrganizer model
%
% Output
% ------
% * VCML file
% * single channel TIF files

% Author: Ivan E. Cao-Berg, Taraz Buck
%
% Copyright (C) 2016-2020 Murphy Lab
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

start_time = tic;
start_cputime = cputime;
disp( mfilename );
disp( 'The estimated running time is about 1 minute. Please wait.' );

% Combine SPHARM framework with vesicle models
combined_models = 'combined_model.mat';
slml2slml({framework_model, vesicle_models}, ...
    struct('output_filename', combined_models, ...
    'selection', [1,1,0;0,0,1]));

options.seed = 639848;
options.targetDirectory = pwd;
options.prefix = 'img';
options.synthesis = 'all';
options.model.spharm_rpdm.synthesis_method = 'random_sampling';
options.model.spharm_rpdm.imageSize = [512, 512, 36];
options.model.spharm_rpdm.synthesis_resolution = [0.10 0.10 0.1];
options.numberOfSynthesizedImages = 1;
options.numberOfGaussianObjects = 5;
options.rendAtStd = 2.0;
options.objstd = options.rendAtStd+0.3;
options.overlapsubsize = 1;
options.overlapthresh = 0;
options.compression = 'lzw';
options.microscope = 'none';
options.sampling.method = 'disc';
options.verbose = true;
options.debug = false;
options.output.tifimages = true;
options.output.indexedimage = true;

% VCML Options
options.output.VCML.writeVCML = true;
options.VCML.downsampling = 1;
options.VCML.translations = {'cell', 'CP'; 'nucleus', 'NU';  ...
    'lamp2_mat_tfr_mat', 'EN'; ...
    'CP_EC', 'PM'; ...
    'CP_EN', 'EM'; ...
    'CP_NU', 'NM'};

options.overlapsubsize = options.VCML.downsampling;
options.overlapthresh = 0;

try
    state = rng( options.seed );
catch err
    rand( 'seed', options.seed ); %#ok<RAND>
end

% Synthesize Image and VCML file
answer = slml2img( { combined_models }, options );
