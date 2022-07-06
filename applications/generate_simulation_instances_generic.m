function answer = generate_simulation_instances_generic( options )
% generate_simulation_instances_generic
%
% Synthesize multiple geometries for simulation and combine with reaction
% network in MCell MDL or Virtual Cell VCML format into multiple simulation
% input files for running using MCell or Virtual Cell.
%
% Input arguments
% ---------------
%
% .. IMPORTANT::
%    Some intermediate and final results are overwritten, and other computations are skipped when the output files already exist. Be aware of this when `output_dir` is set to an existing directory. If you update code or input files, it will often be wise to delete some or all of the contents of this directory prior to running `generate_simulation_instances_generic`.
%
% A structure containing required and optional fields.
% 
% Required fields:
% 
% * `reaction_network_file`: Either a path of a single VCML file or a glob pattern for a single collection of MCell MDL files (default: `'../../data/CBExMinScaled3_20min/CBExMinScaled3_20min.*.mdl'`)
%
% * `output_dir`: The path of a directory into which intermediate and final results will be written
% 
% Optional fields:
% 
% * `framework_model`: The path of a single CellOrganizer framework model (cell and nucleus shapes) (default: `'../../models/3D/spharm/lamp2.mat'`)
% 
% * `vesicle_models`: A cell array of paths of CellOrganizer vesicle models (default: `{'../../models/3D/tfr.mat'}`)
% 
% * `synthesis`: `'all'` to synthesize framework and vesicles or `'framework'` to ignore the vesicle model (default: `'framework'`)
% 
% * `output_image`: Boolean. Whether to additionally produce images of generated geometries (default: `true`)
% 
% * `n_images_to_synthesize`: Number of geometries to sample from `framework_model` to produce distinct simulation inputs (default: `100`)
% 
% * `n_vesicles`: Number of vesicles to sample from each model in `vesicle_models` (default: `[]`)
%     * A value of `[]` samples the number of vesicles from the model's object count distribution
%     * An array of length two containing the minimum and maximum number of objects to be sampled uniformly
% 
% * `vesicle_volume_scale`: Factor by which to scale each vesicle after synthesis (default: `1`)
% 
% * `downsampling`: Downsampling factor to trade disk and computational requirements for geometric precision (default: `1/2`)
% 
% * `vcml_relative_downsampling`: Factor for further downsampling for writing VCML files. VCell simulation data from fewer than 100 simulations at full resolution can fill available disk space on the VCell server, but this is less of a concern when running simulations on your own cluster. (default: `1/2`)
% 
% * `base_seed`: Seed for random number generator (default: `735945`)
% 
% * `n_vesicles_seed_offset`: Offset used to produce seed for a second random number generator for sampling the number of vesicles (default: `34188`)
% 
% * `translations`: N x 2 cell array of strings (first column) to be replaced by other strings (second column) in CellOrganizer-generated VCML, MCell MDL, and SBML (default: `{'cell', 'CP'; 'nuc', 'NU'; 'nucleus', 'NU'; 'lamp2_mat_tfr_mat', 'EN'; 'CP_EC', 'PM'; 'CP_EN', 'EM'; 'CP_NU', 'NM'}`)
% 
% * `simulation_end_time`: Time at which to end simulation in seconds (default: `4000`)
%
% * `simulation_default_time_step`: Initial time step size in seconds for selected simulation method (default: `1e0`)
%
% * `simulation_max_time_step`: Maximum time step size in seconds (default: `4`)
%
% * `simulation_output_time_step`: Output time step size in seconds. Output format varies by simulation method. (default: `100`)
%
% * `simulation_absolute_tolerance`: VCML-specific integrator's absolute error tolerance (default: `1e-8`)
%
% * `simulation_relative_tolerance`: VCML-specific integrator's relative error tolerance (default: `1e-8`)
%
% * `simulation_interaction_radius`: MCell-specific maximum distance for bimolecular reactions to be simulated (default: `0.03`)
%
% * `framework_min_clearance`: double specifying the minimum distance in μm to impose between nucleus and cell after synthesis. -inf to disable. Currently only used for framework meshes assuming corresponding vertices by instance2MCellMDL. (default: `-inf`)
% 
% * `framework_clearance_n_max_filter_rounds`: integer specifying the number of rounds of maximum filter to apply to the projections of cell vertices onto nucleus normals among the immediate neighbors of each vertex. 0 to disable. Currently only used for framework meshes assuming corresponding vertices by instance2MCellMDL. (default: `1`)
% 
% * `intersecting_mesh_object_policy`: string specifying policy for checking framework and objects for intersection and whether to remove objects or reject the synthesized cell entirely. Currently untested for values other than 'ignore'. Currently only used for framework meshes assuming corresponding vertices by instance2MCellMDL. (default: `'ignore'`)
% 
%
% Output
% ------
% % * one valid VCML file with multiple simulations per input VCML model
% % * a set of MCell simulation files per input MCell model

% Developer notes
% ---------------
% 
% * Only certain options for `slml2img` are copied from `options` to `s2i_options`. Add options for which support is desired when using this utility around the lines `% Default arguments` and `% Copy options fields for use by slml2img`.

% 
% * ``:  (default: `''`)

% Author: Taraz Buck, Ivan E. Cao-Berg
%
% Copyright (C) 2016-2022 Murphy Lab
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

% Default arguments
options_defaults = struct();
options_defaults.framework_model = '../../models/3D/spharm/lamp2.mat';
options_defaults.vesicle_models = {'../../models/3D/tfr.mat'};
options_defaults.synthesis = 'framework';
%{
if strcmp(options_defaults.synthesis, 'all')
    % options_defaults.reaction_network_file = '../../data/CBExMinScaled3EN_20min.vcml';
    options_defaults.reaction_network_file = '../../data/CBExMinScaled3EN20min/CBExMinScaled3EN20min.*.mdl';
else
    % options_defaults.reaction_network_file = '../../data/CBExMinScaled3_20min.vcml';
    options_defaults.reaction_network_file = '../../data/CBExMinScaled3_20min/CBExMinScaled3_20min.*.mdl';
end
%}
% options_defaults.output_vcml = true;
% options_defaults.output_mcell = false;
options_defaults.output_image = true;

options_defaults.n_images_to_synthesize = 100;
options_defaults.n_vesicles = [];
options_defaults.vesicle_volume_scale = 1;

options_defaults.downsampling = 1/2;
options_defaults.vcml_relative_downsampling = 1/2;

options_defaults.base_seed = 735945;
options_defaults.n_vesicles_seed_offset = 34188;

options_defaults.translations = {'cell', 'CP'; 'nuc', 'NU'; 'nucleus', 'NU'; 'lamp2_mat_tfr_mat', 'EN'; 'CP_EC', 'PM'; 'CP_EN', 'EM'; 'CP_NU', 'NM'};

options_defaults.use_old_min_naming = false;

options_defaults.simulation_end_time = 4000;
options_defaults.simulation_default_time_step = 1e0;
options_defaults.simulation_max_time_step = 4;
options_defaults.simulation_output_time_step = 100;
% Virtual Cell-specific
options_defaults.simulation_absolute_tolerance = 1e-8;
options_defaults.simulation_relative_tolerance = 1e-8;
% MCell-specific
options_defaults.simulation_interaction_radius = 0.03;

options_defaults.framework_min_clearance = -inf;
options_defaults.framework_clearance_n_max_filter_rounds = 1;
options_defaults.intersecting_mesh_object_policy_default = 'ignore';

options_required = struct();
options_required.reaction_network_file = [];
options_required.output_dir = [];

options = ml_initparam( options, options_defaults, options_required );

% (framework_model|vesicle_models|synthesis|reaction_network_file|output_image|n_images_to_synthesize|n_vesicles|vesicle_volume_scale|downsampling|vcml_relative_downsampling|base_seed|translations|output_dir)

%{
% The framework (cell and nucleus shapes) will be synthesized using framework_model
framework_model = options.framework_model;
vesicle_models = options.vesicle_models;
% Set to 'all' to enable object generation
synthesis = options.synthesis;
reaction_network_file = options.reaction_network_file;
output_formats = {};
if options.output_vcml; output_formats{end+1} = 'vcml'; end;
if options.output_mcell; output_formats{end+1} = 'mcell'; end;
% if options.output_sbml; output_formats{end+1} = 'sbml'; end;
if options.output_image; output_formats{end+1} = 'image'; end;
base_seed = options.base_seed;
n_images_to_synthesize = options.n_images_to_synthesize;
n_vesicles = options.n_vesicles;
vesicle_volume_scale = options.vesicle_volume_scale;
downsampling = options.downsampling;
vcml_relative_downsampling = options.vcml_relative_downsampling;
%}
fprintf('vcml_relative_downsampling = %f\n', options.vcml_relative_downsampling);

% error('Not implemented');
if ~isdeployed()
    %{
    current_path = which(mfilename);
    [current_path, filename, extension] = fileparts( current_path );
    current_path = [current_path, filesep, mfilename, filesep];
    %}
    current_path = options.output_dir;
    current_path
    mkdir_recursive(current_path);
    cd(current_path);
end

start_time = tic;
start_cputime = cputime;
disp( mfilename );
% disp( ['The estimated running time is about ', num2str(options.n_images_to_synthesize * 20/16, '%.1f'), ' minutes. Please wait.'] );
% disp( ['The estimated running time is about ', num2str(options.n_images_to_synthesize * 196/60, '%.1f'), ' minutes. Please wait.'] );
disp( ['The estimated running time is about ', num2str(options.n_images_to_synthesize * 24478.721/100/60, '%.1f'), ' minutes. Please wait.'] );


debug_single_geometry = false;
% debug_single_geometry = true
debug_quick_geometry = false;
% debug_quick_geometry = true
debug_always_overwrite = false;
% debug_always_overwrite = true

debug_quick_simulation = false;
% debug_quick_simulation = true
debug_downsampling = false;
% debug_downsampling = true

debug_output_count_offset = 0;
% debug_output_count_offset = -2

debug_zero_objects = false;
% debug_zero_objects = true


% Combine SPHARM framework with vesicle models
combined_models = 'combined_model';
chunk_filename = combined_models;
[can_start, final_name, final_exists, temp_name] = chunk_start(chunk_filename, '.mat', true, 5 * 60);
combined_models = final_name;
if debug_always_overwrite
    can_start = true;
end
if can_start
    slml2slml({options.framework_model, options.vesicle_models}, struct('output_filename', chunk_filename, 'selection', [1,1,0;0,0,1]));
    chunk_finish(chunk_filename);
end
if ~exist(final_name, 'file')
    error('"%s" does not exist', final_name);
end

if options.vesicle_volume_scale ~= 1 && strcmp(options.synthesis, 'all')
    % Scale vesicles so fewer have to be generated to occupy space, reducing synthesis time and mesh file size
    error('Vesicle scaling not working');
    [temp2, temp3, temp4] = fileparts(combined_models);
    temp3 = [temp3, sprintf('_vesvolscale%.3f', options.vesicle_volume_scale)];
    % combined_models2 = fullfile(temp2, strcat(temp3, temp4));
    combined_models2 = fullfile(temp2, temp3);
    chunk_filename = combined_models2;
    [can_start, final_name, final_exists, temp_name] = chunk_start(chunk_filename, '.mat', true, 5 * 60);
    combined_models2 = final_name;
    if debug_always_overwrite
        can_start = true;
    end
    if can_start
        temp = load(combined_models);
        temp = temp.model;
        proteinModel = temp.proteinModel
        % dbstack, keyboard % Debug
        vesicle_volume_scale_cube_root = options.vesicle_volume_scale^(1/3);
        % http://repositories.compbio.cs.cmu.edu/cellorganizer/cellorganizer3/blob/dev/utilities/3D/vesicles/3D/gmm_objsizefit_percell.m
        % Log-normal mean and standard deviation
        proteinModel.size.x.mu = proteinModel.size.x.mu + log(vesicle_volume_scale_cube_root);
        proteinModel.size.x.sigma = proteinModel.size.x.sigma + log(vesicle_volume_scale_cube_root);
        % Conditional normal
        % Mean of y size
        proteinModel.size.y_x.a1 = proteinModel.size.y_x.a1 + log(vesicle_volume_scale_cube_root);
        % Standard deviation of y size
        proteinModel.size.y_x.a2 = proteinModel.size.y_x.a2 + log(vesicle_volume_scale_cube_root);
        % Unknown, appears to be unaffected by x and y scaling
        proteinModel.size.y_x.b1 = proteinModel.size.y_x.b1 + 0;
        proteinModel.size.y_x.b2 = proteinModel.size.y_x.b2 + 0;
        % Conditional normal
        % Mean of z size
        proteinModel.size.z_x.a1 = proteinModel.size.z_x.a1 + log(vesicle_volume_scale_cube_root);
        % Standard deviation of z size
        proteinModel.size.z_x.a2 = proteinModel.size.z_x.a2 + log(vesicle_volume_scale_cube_root);
        % Unknown, appears to be unaffected by x and z scaling
        proteinModel.size.z_x.b1 = proteinModel.size.z_x.b1 + 0;
        proteinModel.size.z_x.b2 = proteinModel.size.z_x.b2 + 0;
        temp.proteinModel = proteinModel;
        temp = struct('model', temp);
        save(combined_models2, '-struct', 'temp');
        
        chunk_finish(chunk_filename);
    end
    combined_models = combined_models2;
end

n_images_to_synthesize_for_debugging = options.n_images_to_synthesize;
n_images_to_skip_for_debugging = 0;

mcell_extensions = {'.mdl', '.mcell'};
vcml_extensions = {'.vcml'};
mcell_input_extensions = mcell_extensions;
% mcell_input_extensions = [mcell_input_extensions, {'.net'}];
vcml_input_extensions = vcml_extensions;
% vcml_input_extensions = [vcml_input_extensions, {'', '.net'}];
reaction_network_input_extensions = [mcell_input_extensions, vcml_input_extensions];

[reaction_network_file_path, reaction_network_file_name, reaction_network_file_ext] = fileparts(options.reaction_network_file);

reaction_network_file_name2 = reaction_network_file_name;
if any(strcmpi(reaction_network_file_ext, mcell_extensions))
    reaction_network_file_name2 = strrep(reaction_network_file_name2, '.*', '');
end
prefix = 'img';
if any(strcmpi(reaction_network_file_ext, reaction_network_input_extensions))
    prefix = sprintf('%s_%s', prefix, reaction_network_file_name2);
end

% write_SBMLSpatial = any(strcmp(output_formats, 'sbml')) && (isempty(options.reaction_network_file) || any(strcmpi(reaction_network_file_ext, {'.xml'})));
% write_VCML = any(strcmp(output_formats, 'vcml')) && (isempty(options.reaction_network_file) || any(strcmpi(reaction_network_file_ext, vcml_input_extensions)));
% write_MCellMDL = any(strcmp(output_formats, 'mcell')) && (isempty(options.reaction_network_file) || any(strcmpi(reaction_network_file_ext, mcell_input_extensions)));
% write_images = any(strcmp(output_formats, 'image'));
write_SBMLSpatial = any(strcmpi(reaction_network_file_ext, {'.xml'}));
write_VCML = any(strcmpi(reaction_network_file_ext, vcml_input_extensions));
write_MCellMDL = any(strcmpi(reaction_network_file_ext, mcell_input_extensions));
write_images = options.output_image;

if debug_downsampling
    n_images_to_synthesize_for_debugging = 1;
    % n_images_to_synthesize_for_debugging = 2;
    
    % options.downsampling = 1/8;
    options.downsampling = 1;
end

n_vesicles_mode = 'model';
n_vesicles_min = [];
n_vesicles_max = [];
if isempty(options.n_vesicles)
    % Default
elseif isnumeric(options.n_vesicles) && length(options.n_vesicles) == 2
    n_vesicles_mode = 'uniform';
    n_vesicles_min = options.n_vesicles(1);
    n_vesicles_max = options.n_vesicles(2);
else
    error('`options.n_vesicles` must be `[]` or an array of two integers');
end

%{
% % Unused for this example, but set to enable object generation
% Was `n_gaussian_objects_min`
n_vesicles_min = 1;
% n_vesicles_max = 200;
% n_vesicles_max = 1500;
% n_vesicles_max = 3000;
% Attempt to fill between 0 and cytoplasm_max_mean_fill_fraction of cytoplasm. Intersecting objects will be discarded.
cytoplasm_max_mean_fill_fraction = 0.5;
% Computed from data from these 102 samples ('generate_simulation_instances_min.20201006').
cytoplasm_volume_mean = 975.6973273003645;
vesicle_volume_mean = 0.035302936430848024 * vesicle_volume_scale;
n_vesicles_max = ceil(cytoplasm_max_mean_fill_fraction * cytoplasm_volume_mean / vesicle_volume_mean);
n_vesicles_min = ceil(n_vesicles_max / 10);
% fprintf('n_vesicles_max = %f\n', n_vesicles_max);
% fprintf('n_vesicles_min = %f\n', n_vesicles_min);
%}
if debug_zero_objects
    n_vesicles_mode = 'uniform'
    n_vesicles_min = 0
    n_vesicles_max = 200
    options.base_seed = 735945;
end
if debug_quick_geometry
    n_vesicles_mode = 'uniform'
    n_vesicles_min = 25
    n_vesicles_max = 25
end
function [result] = n_vesicles_function(given_seed)
    switch n_vesicles_mode
        case 'model'
            result = [];
        case 'uniform'
            if n_vesicles_min == n_vesicles_max
                result = n_vesicles_min;
            else
                seed2 = given_seed + options.n_vesicles_seed_offset;
                try
                    state2 = rng( seed2 );
                catch err
                    state2 = rand( 'seed', seed2 ); %#ok<RAND>
                end
                result = randi([n_vesicles_min, n_vesicles_max]);
            end
            try
                rng( state2 );
            catch err
                rand( 'state', state2 ); %#ok<RAND>
            end
    end
end

if debug_single_geometry || debug_quick_geometry
    n_images_to_synthesize_for_debugging = 1;
    % n_images_to_skip_for_debugging = 0;
    n_images_to_skip_for_debugging = 2;
end
if debug_zero_objects
    n_images_to_synthesize_for_debugging = 1;
    n_images_to_skip_for_debugging = 0;
end

if debug_quick_geometry
    % n_images_to_synthesize_for_debugging = 2;
    
    options.downsampling = 1/2; % Debug
end

function [result_seed_prefix, result_seed_temporary_results, result_chunk_filename, result_can_start, result_final_name, result_final_exists, result_temp_name] = custom_chunk_start(given_seed, given_prefix, given_suffix, given_extension)
    if nargin < 4
        given_extension = '.mat';
    end
    result_seed_prefix = sprintf('%s_%010d%s', given_prefix, given_seed, given_suffix);
    result_seed_temporary_results = [pwd, filesep, result_seed_prefix, '_temp'];
    result_chunk_filename = result_seed_prefix;
    
    [result_can_start, result_final_name, result_final_exists, result_temp_name] = chunk_start(result_chunk_filename, given_extension);
    
    if debug_always_overwrite
        result_can_start = true;
    end
end


%{
function [result] = mdl_count_function(given_prefix, given_suffix)
    mdl_pattern = [given_prefix, '*', given_suffix, filesep, 'cell1', filesep, 'cell.geometry.mdl'];
    result = length(dir(mdl_pattern)) + debug_output_count_offset;
end

function [result] = vcml_count_function(given_prefix, given_suffix)
    vcml_pattern = [given_prefix, '*', given_suffix, filesep, 'cell1', filesep, 'cell.vcml'];
    result = length(dir(vcml_pattern)) + debug_output_count_offset;
end
%}

%{
function [result] = temp_count_function(given_prefix, given_suffix)
    result = 0;
    temp_pattern = [given_prefix, '*', given_suffix, '.tmp'];
    result = result + length(dir(temp_pattern));
end
%}

function [result] = generic_output_count_function(given_prefix, given_suffix, given_suffix2, given_extensions)
    % given_prefix, given_suffix, given_suffix2, given_extensions % Debug
    result = 0;
    for i=1:length(given_extensions)
        dir_pattern = [given_prefix, '*', given_suffix];
        given_suffix2_with_extension = [given_suffix2, given_extensions{i}];
        combined_pattern = [dir_pattern, filesep, given_suffix2_with_extension];
        % tmp_pattern = [dir_pattern, '.tmp'];
        % result = result + length(combined_pattern_matches);
        combined_pattern_matches = dir(combined_pattern);
        % dir_pattern, given_suffix2_with_extension, combined_pattern, combined_pattern_matches % Debug
        for j=1:length(combined_pattern_matches)
            tmp_filename = [combined_pattern(1:end - (length(given_suffix2_with_extension) + 1)), '.tmp'];
            if ~exist(tmp_filename, 'file')
                result = result + 1;
            end
        end
        %{
        result % Debug
        % error('Check for .tmp file for each entry!');
        if any(strcmp(given_extensions, '.vcml')) % Debug
            error('Not implemented')
        end
        %}
    end
    result = result + debug_output_count_offset;
end

function [result] = mcell_count_function(given_prefix, given_suffix)
    result = generic_output_count_function(given_prefix, given_suffix, ['cell1', filesep, 'cell.geometry'], mcell_extensions);
end

function [result] = vcml_count_function(given_prefix, given_suffix)
    result = generic_output_count_function(given_prefix, given_suffix, ['cell1', filesep, 'cell'], vcml_extensions);
end

function [result] = output_count_function(given_prefix, given_suffix)
    mcell_count = mcell_count_function(given_prefix, given_suffix);
    vcml_count = vcml_count_function(given_prefix, given_suffix);
    result = [];
    if write_MCellMDL
        result(end+1) = mcell_count;
    end
    if write_VCML
        result(end+1) = vcml_count;
    end
    result = min(result);
end


suffix = '';
if debug_quick_simulation
    options.simulation_end_time = options_defaults.simulation_end_time / 10;
end
if options.use_old_min_naming
    suffix = sprintf('%s_%04dsec', suffix, options.simulation_end_time);
else
    suffix = sprintf('%s_%.4esec', suffix, options.simulation_end_time);
end

s2i_options = struct();
s2i_options.clean_synthetic_instances = false;
s2i_options.seed = options.base_seed;
s2i_options.targetDirectory = pwd;
s2i_options.synthesis = options.synthesis;
s2i_options.model.spharm_rpdm.synthesis_method = 'random_sampling';
s2i_options.model.spharm_rpdm.imageSize = [205, 205, 18];
s2i_options.numberOfSynthesizedImages = 1;
s2i_options.output.remove_mesh_intersections = true;

if strcmp(s2i_options.synthesis, 'all')
    % I think this is done because slml2img or a function it calls resamples the framework by a factor of 4 when synthesizing protein distributions
    downsampling2 = options.downsampling * 1/4;
elseif strcmp(s2i_options.synthesis, 'framework')
    downsampling2 = options.downsampling * 1;
end
if write_VCML
    downsampling2 = downsampling2 * options.vcml_relative_downsampling;
end
fprintf('downsampling2 = %f\n', downsampling2);

% s2i_options.numberOfGaussianObjects = options.n_vesicles;

s2i_options.rendAtStd = 1;
% s2i_options.rendAtStd = 2;
% s2i_options.objstd = s2i_options.rendAtStd+0.3;
s2i_options.overlapsubsize = 1;
s2i_options.overlapthresh = 0;
s2i_options.oobbuffer = 0.1;
if debug_quick_geometry
    % s2i_options.oobbuffer = 1; % Debug
end
s2i_options.compression = 'lzw';
s2i_options.microscope = 'none';
s2i_options.sampling.method = 'disc';
s2i_options.verbose = true;
s2i_options.debug = false;
s2i_options.output.tifimages = write_images;
s2i_options.output.shape_space_coords = true;
s2i_options.output.OMETIFF = write_images;
s2i_options.output.indexedimage = write_images;
s2i_options.output.meshes = true;
s2i_options.SBML.downsampling = downsampling2;
s2i_options.SBML.spatial_use_compression = true;
s2i_options.SBML.spatial_use_analytic_meshes = true;
s2i_options.SBML.spatial_image = false;
s2i_options.SBML.spatial_vcell_compatible = false;
s2i_options.SBML.include_ec = true;
s2i_options.SBML.ec_scale = 1.25;

s2i_options.output.SBMLSpatial = write_SBMLSpatial;
s2i_options.output.VCML.writeVCML = write_VCML;
s2i_options.output.writeMCellMDL = write_MCellMDL;

s2i_options.VCML.downsampling = downsampling2;

if isempty(options.reaction_network_file)
%{
elseif strcmpi(reaction_network_file_ext, '.net')
    s2i_options.NET.filename = options.reaction_network_file;
    s2i_options.NET.units.length = 'm';
    s2i_options.NET.units.concentration = 'mol.m-3';
    s2i_options.NET.useImageAdjacency = false;
%}
elseif write_VCML
    s2i_options.VCML.input_filename = options.reaction_network_file;
    s2i_options.NET.useImageAdjacency = false;
elseif write_MCellMDL
    s2i_options.MCellMDL.input_filename_pattern = options.reaction_network_file;
    s2i_options.NET.useImageAdjacency = false;
else
    error('reaction network file of type ''%s'' not supported', reaction_network_file_ext);
end



if options.use_old_min_naming
    prefix = sprintf('%s_%2.2f', prefix, downsampling2);
else
    prefix = sprintf('%s_%.2e', prefix, downsampling2);
end
s2i_options.prefix = prefix;

s2i_options.NET.translations = options.translations;

s2i_options.VCML.translations = options.translations;
s2i_options.VCML.end_time = options.simulation_end_time;
s2i_options.VCML.default_time_step = options.simulation_default_time_step;
s2i_options.VCML.max_time_step = options.simulation_max_time_step;
s2i_options.VCML.output_time_step = options.simulation_output_time_step;
s2i_options.VCML.absolute_tolerance = options.simulation_absolute_tolerance;
s2i_options.VCML.relative_tolerance = options.simulation_relative_tolerance;

s2i_options.SBML.translations = options.translations;
s2i_options.SBML.end_time = options.simulation_end_time;
s2i_options.SBML.default_time_step = options.simulation_default_time_step;
s2i_options.SBML.max_time_step = options.simulation_max_time_step;
s2i_options.SBML.output_time_step = options.simulation_output_time_step;

s2i_options.MCellMDL.translations = options.translations;
s2i_options.MCellMDL.end_time = options.simulation_end_time;
s2i_options.MCellMDL.default_time_step = options.simulation_default_time_step;
s2i_options.MCellMDL.max_time_step = options.simulation_max_time_step;
s2i_options.MCellMDL.output_time_step = options.simulation_output_time_step;
s2i_options.MCellMDL.interaction_radius = options.simulation_interaction_radius;

% Copy options fields for use by slml2img

s2i_options.framework_min_clearance = options.framework_min_clearance;
s2i_options.framework_clearance_n_max_filter_rounds = options.framework_clearance_n_max_filter_rounds;
s2i_options.intersecting_mesh_object_policy_default = options.intersecting_mesh_object_policy_default;





% Produce files for at least n_images_to_synthesize_for_debugging simulations
seed_base = options.base_seed + n_images_to_skip_for_debugging;

seed_index = 0;
seeds_finished = containers.Map('KeyType', 'uint32', 'ValueType', 'logical');
fprintf('\nOutput count: %i\n\n', output_count_function(prefix, suffix));
% while output_count_function(prefix, suffix) < n_images_to_synthesize_for_debugging
while (output_count_function(prefix, suffix) < n_images_to_synthesize_for_debugging && ~debug_zero_objects) || (debug_zero_objects && seed_index == 0)
    seed = seed_base + seed_index;
    seed_index = seed_index + 1;
    
    seeds_finished(seed) = false;
    [seed_prefix, seed_temporary_results, chunk_filename, can_start, final_name, final_exists, temp_name] = custom_chunk_start(seed, prefix, suffix);
    
    if debug_zero_objects
        can_start = true;
    end
    
    if ~can_start
        if final_exists
            seeds_finished(seed) = true;
        end
        continue;
    end

    slml2img_start_time = tic;
    slml2img_start_cputime = cputime;
    
    s2i_options.seed = seed;
    s2i_options.prefix = seed_prefix;
    s2i_options.temporary_results = seed_temporary_results;
    
    s2i_options.numberOfGaussianObjects = n_vesicles_function(seed);
    fprintf('\noptions.numberOfGaussianObjects: %i\n\n', s2i_options.numberOfGaussianObjects);
    
    try
        state = rng( s2i_options.seed );
    catch err
        rand( 'seed', s2i_options.seed ); %#ok<RAND>
    end
    answer = slml2img( {combined_models}, s2i_options );

    slml2img_elapsed_time = toc(slml2img_start_time);
    slml2img_elapsed_cputime = cputime - slml2img_start_cputime;
    fprintf('\nslml2img for %s took %.3f s (%.3f s CPU time)\n\n', seed_prefix, slml2img_elapsed_time, slml2img_elapsed_cputime);
    
    if answer
        empty = [];
        save(final_name, 'empty');
        seeds_finished(seed) = true;
    end
    
    chunk_finish(chunk_filename);
    fprintf('\nOutput count: %i\n\n', output_count_function(prefix, suffix));
end


seeds_finished_values = all(cell2mat(seeds_finished.values()));
if all(seeds_finished_values) && write_VCML
    % Combine all VCML files into one for convenience
    
    combined_vcml_name = sprintf('%s_%010d_%04d%s', prefix, options.base_seed, n_images_to_synthesize_for_debugging, suffix)
    [combined_can_start, combined_final_name, combined_final_exists, combined_temp_name] = chunk_start(combined_vcml_name, '.vcml');
    if debug_always_overwrite
        combined_can_start = true;
    end
    
    if combined_can_start
        fprintf('\nCombining VCML files\n');
        combining_start_time = tic;
        combining_start_cputime = cputime;
        
        all_vcml = [];
        all_vcml_biomodel = [];
        for seed = seeds_finished_values
            [~, ~, seed_chunk_filename, ~, ~, ~, ~] = custom_chunk_start(seed, prefix, suffix);
            
            vcml_name = [seed_chunk_filename, filesep, 'cell1', filesep, 'cell.vcml'];
            seed_vcml = xmlread(vcml_name);
            if isempty(all_vcml)
                all_vcml = seed_vcml;
                all_vcml_biomodel = all_vcml.getElementsByTagName('BioModel');
                if all_vcml_biomodel.getLength() ~= 1
                    error('VCML file should have one BioModel');
                end
                all_vcml_biomodel = all_vcml_biomodel.item(0);
            else
                seed_vcml_biomodel = seed_vcml.getElementsByTagName('BioModel');
                if seed_vcml_biomodel.getLength() ~= 1
                    error('VCML file should have one BioModel');
                end
                seed_vcml_biomodel = seed_vcml_biomodel.item(0);
                seed_vcml_simulationspecs = seed_vcml_biomodel.getElementsByTagName('SimulationSpec');
                if seed_vcml_simulationspecs.getLength() == 0
                    error('VCML file should have at least one SimulationSpec');
                end
                for k = 1:seed_vcml_simulationspecs.getLength()
                    seed_vcml_simulationspec = seed_vcml_simulationspecs.item(k-1);
                    all_vcml_biomodel.appendChild(all_vcml.importNode(seed_vcml_simulationspec, true));
                end
            end
        end
        
        xmlwrite(combined_final_name, all_vcml);
        
        chunk_finish(combined_vcml_name);
        
        combining_elapsed_time = toc(combining_start_time);
        combining_elapsed_cputime = cputime - combining_start_cputime;
        fprintf('\nCombining VCML files took %.3f s (%.3f s CPU time)\n\n', combining_elapsed_time, combining_elapsed_cputime);
    end
end



elapsed_time = toc(start_time);
elapsed_cputime = cputime - start_cputime;
fprintf('\n%s took %.3f s (%.3f s CPU time)\n\n', mfilename, elapsed_time, elapsed_cputime);

end%generate_simulation_instances_SarmaGhosh2012
