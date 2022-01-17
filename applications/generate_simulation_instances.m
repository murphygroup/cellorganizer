function answer = generate_simulation_instances( options )
% generate_simulation_instances
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
framework_model = '../../models/3D/spharm/lamp2.mat';
vesicle_models = '../../models/3D/tfr.mat';
% net_file = '../../data/SavageEtAl2012ReactionDiffusion.net';
net_file = '../../data/LotkaVolterra.net';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT MODIFY THIS BLOCK
if ~isdeployed()
  current_path = which(mfilename);
  [current_path, filename, extension] = fileparts( current_path );
  % save_path = [pwd, filesep, mfilename, filesep];
  % current_path = [pwd, filesep, mfilename, filesep];
  current_path = [current_path, filesep, mfilename, filesep];
  mkdir_recursive(current_path);
  cd(current_path);
end

start_time = tic;
start_cputime = cputime;
disp( mfilename );
disp( 'The estimated running time is about 3.5 hours. Please wait.' );

% Combine SPHARM framework with vesicle models
combined_models = 'combined_model.mat';
slml2slml({framework_model, vesicle_models}, struct('output_filename', combined_models, 'selection', [1,1,0;0,0,1]));

base_seed = 315200;
n_images_to_synthesize = 100;
[net_file_path, net_file_name, net_file_ext] = fileparts(net_file);
prefix = ['img_', net_file_name];

options.seed = base_seed;
options.targetDirectory = pwd;
% options.targetDirectory = save_path;
options.prefix = prefix;
options.synthesis = 'all';
options.model.spharm_rpdm.synthesis_method = 'random_sampling';
options.model.spharm_rpdm.imageSize = [205, 205, 18];
options.numberOfSynthesizedImages = 1;
% options.numberOfGaussianObjects = 10;
options.numberOfGaussianObjects = 50;
% options.numberOfGaussianObjects = 100;
options.rendAtStd = 1;
options.overlapthresh = 1;
options.compression = 'lzw';
options.microscope = 'none';
options.sampling.method = 'disc';
options.verbose = true;
options.debug = false;
options.output.tifimages = true;
options.output.shape_space_coords = true;
options.output.OMETIFF = true;
options.output.indexedimage = true;
% options.output.SBMLSpatial = true;
% options.SBML.downsampling = [1/4, 1/4, 1];
options.SBML.spatial_use_compression = true;
options.SBML.spatial_use_analytic_meshes = true;
options.SBML.spatial_image = false;
options.SBML.spatial_vcell_compatible = false;
options.output.SBMLSpatial = true;
% options.output.writeVCML = true;
options.VCML.downsampling = 14 / 73;
options.output.writeMCellMDL = true;
options.oobbuffer = 0.1;

options.NET.filename = net_file;
options.NET.units.length = 'm';
options.NET.units.concentration = 'mol.m-3';
options.NET.translations = {'cell', 'CP'; 'nuc', 'NU'; 'nucleus', 'NU'; 'lamp2_mat_tfr_mat', 'EN'; 'CP_EC', 'PM'; 'CP_EN', 'EM'; 'CP_NU', 'NM'};
options.NET.useImageAdjacency = false;

%{
options.VCML.end_time = 36;
options.VCML.default_time_step = 0.05;
options.VCML.max_time_step = 1;
%}
%{
options.VCML.end_time = 360;
options.VCML.default_time_step = 0.05;
options.VCML.max_time_step = 1;
%}
options.VCML.default_time_step = 3.3333333333333333e-7;
options.VCML.end_time = options.VCML.default_time_step * 6e4;
options.VCML.max_time_step = options.VCML.default_time_step * 100;
options.VCML.output_time_step = options.VCML.default_time_step * 120;

options.MCellMDL.end_time = options.VCML.end_time;
options.MCellMDL.default_time_step = options.VCML.default_time_step;
options.MCellMDL.max_time_step = options.VCML.max_time_step;
options.MCellMDL.output_time_step = options.VCML.output_time_step;

options.MCellMDL.interaction_radius = 0.03;

options.MCellMDL.use_reaction_rate_hack = true;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


all_seeds = base_seed + (0:n_images_to_synthesize - 1);
for seed = all_seeds
    seed_prefix = sprintf('%s_%010d', prefix, seed);
    seed_temporary_results = [pwd, filesep, seed_prefix, '_temp'];
    chunk_filename = seed_prefix;
    
    [can_start, final_name, final_exists, temp_name] = chunk_start(chunk_filename, '.mat');
    
    if ~can_start
        continue;
    end

    slml2img_start_time = tic;
    slml2img_start_cputime = cputime;
    
    options.seed = seed;
    options.prefix = seed_prefix;
    options.temporary_results = seed_temporary_results;
    
    try
        state = rng( options.seed );
    catch err
        rand( 'seed', options.seed ); %#ok<RAND>
    end
    answer = slml2img( {combined_models}, options );
    
    empty = [];
    save(final_name, 'empty');
    
    chunk_finish(chunk_filename);

    slml2img_elapsed_time = toc(slml2img_start_time);
    slml2img_elapsed_cputime = cputime - slml2img_start_cputime;
    fprintf('\nslml2img for %s took %.3f s (%.3f s CPU time)\n\n', seed_prefix, slml2img_elapsed_time, slml2img_elapsed_cputime);
end

elapsed_time = toc(start_time);
elapsed_cputime = cputime - start_cputime;
fprintf('\n%s took %.3f s (%.3f s CPU time)\n\n', mfilename, elapsed_time, elapsed_cputime);

end%generate_simulation_instances
