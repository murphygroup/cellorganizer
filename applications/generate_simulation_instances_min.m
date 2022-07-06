function answer = generate_simulation_instances_min( options )
% generate_simulation_instances_min
%
% Synthesize multiple geometries for simulation and combine with reaction
% network in MCell MDL format into directories ready for MCell simulation.
%
% Input
% -----
% * a valid CellOrganizer model
% * one or more valid VCML reaction network models with Virtual Cell's 
%   default units
%
% Output
% ------
% * one or more valid VCML files with multiple simulations

% Author: Taraz Buck, Ivan E. Cao-Berg
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
framework_model = '../../models/3D/spharm/lamp2.mat';
vesicle_models = '../../models/3D/tfr.mat';
% Set to 'all' to enable object generation
% synthesis = 'framework';
synthesis = 'all';

% reaction_network_files is a cell array of cell arrays. Each entry in the first level has a different random seed, so VCML files given in a second level will share geometries. This allows reactions to be compared by controlling for geometric differences.
reaction_network_files = {};
reaction_network_files{end+1} = {};

reaction_network_files{end}{end+1} = '../../data/CBExMinScaled3EN20min/CBExMinScaled3EN20min.*.mdl';

output_formats = {};
% output_formats{end+1} = 'vcml';
output_formats{end+1} = 'mcell'
% output_formats{end+1} = 'sbml'
output_formats{end+1} = 'image';

base_seed = 735860;
n_images_to_synthesize = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT MODIFY THIS BLOCK
if ~isdeployed()
  current_path = which(mfilename);
  [current_path, filename, extension] = fileparts( current_path );
  current_path = [current_path, filesep, mfilename, filesep];
  current_path
  mkdir_recursive(current_path);
  cd(current_path);
end

start_time = tic;
start_cputime = cputime;
disp( mfilename );
% disp( ['The estimated running time is about ', num2str(n_images_to_synthesize * 20/16, '%.1f'), ' minutes. Please wait.'] );
% disp( ['The estimated running time is about ', num2str(n_images_to_synthesize * 196/60, '%.1f'), ' minutes. Please wait.'] );
disp( ['The estimated running time is about ', num2str(n_images_to_synthesize * 24478.721/100/60, '%.1f'), ' minutes. Please wait.'] );

% Combine SPHARM framework with vesicle models
combined_models = 'combined_model.mat';
chunk_filename = combined_models;
[can_start, final_name, final_exists, temp_name] = chunk_start(chunk_filename, '.mat', true, 5 * 60);
if can_start
    slml2slml({framework_model, vesicle_models}, struct('output_filename', chunk_filename, 'selection', [1,1,0;0,0,1]));
    chunk_finish(chunk_filename);
end

n_images_to_synthesize = n_images_to_synthesize;

% VCell simulation data from fewer than 100 simulations can fill available disk space on the VCell server
downsampling = 1/2;

% % Unused for this example, but set to enable object generation
% n_gaussian_objects_function = @(given_seed)50;
n_gaussian_objects_max = 200;
function [result] = n_gaussian_objects_function(given_seed)
    seed2 = given_seed + 34188;
    try
        state2 = rng( seed2 );
    catch err
        state2 = rand( 'seed', seed2 ); %#ok<RAND>
    end
    result = randi([0, n_gaussian_objects_max]);
    try
        rng( state2 );
    catch err
        rand( 'state', state2 ); %#ok<RAND>
    end
end

function [result_seed_prefix, result_seed_temporary_results, result_chunk_filename, result_can_start, result_final_name, result_final_exists, result_temp_name] = custom_chunk_start(given_seed, given_prefix, given_suffix, given_extension)
    if nargin < 4
        given_extension = '.mat';
    end
    result_seed_prefix = sprintf('%s_%010d%s', given_prefix, given_seed, given_suffix);
    result_seed_temporary_results = [pwd, filesep, result_seed_prefix, '_temp'];
    result_chunk_filename = result_seed_prefix;
    
    [result_can_start, result_final_name, result_final_exists, result_temp_name] = chunk_start(result_chunk_filename, given_extension);
end


for i = 1:length(reaction_network_files)
    base_seed_reaction_network_files = reaction_network_files{i};
    for j = 1:length(base_seed_reaction_network_files)
        reaction_network_file = base_seed_reaction_network_files{j};
        
        % Base seed changes with i, so group reaction network files to use the same geometries
        reaction_network_file_base_seed = base_seed + (i - 1) * n_images_to_synthesize;
        
        [reaction_network_file_path, reaction_network_file_name, reaction_network_file_ext] = fileparts(reaction_network_file);
        suffix = '';

        options = struct();
        options.clean_synthetic_instances = false;
        options.seed = reaction_network_file_base_seed;
        options.targetDirectory = pwd;
        options.synthesis = synthesis;
        options.model.spharm_rpdm.synthesis_method = 'random_sampling';
        options.model.spharm_rpdm.imageSize = [205, 205, 18];
        options.numberOfSynthesizedImages = 1;
        options.output.remove_mesh_intersections = true;

        if strcmp(options.synthesis, 'all')
            downsampling = downsampling * 1/4;
        elseif strcmp(options.synthesis, 'framework')
            downsampling = downsampling * 1;
        end
        fprintf('downsampling = %f\n', downsampling);
        
        % options.numberOfGaussianObjects = n_gaussian_objects;

        options.rendAtStd = 1;
        % options.rendAtStd = 2;
        % options.objstd = options.rendAtStd+0.3;
        options.overlapsubsize = 1;
        options.overlapthresh = 0;
        options.oobbuffer = 0.1;
        options.compression = 'lzw';
        options.microscope = 'none';
        options.sampling.method = 'disc';
        options.verbose = true;
        options.debug = false;
        options.output.tifimages = any(strcmp(output_formats, 'image'));
        options.output.shape_space_coords = true;
        options.output.OMETIFF = any(strcmp(output_formats, 'image'));
        options.output.indexedimage = any(strcmp(output_formats, 'image'));
        options.SBML.downsampling = downsampling;
        options.SBML.spatial_use_compression = true;
        options.SBML.spatial_use_analytic_meshes = true;
        options.SBML.spatial_image = false;
        options.SBML.spatial_vcell_compatible = false;
        options.SBML.include_ec = true;
        options.SBML.ec_scale = 1.25;

        options.output.SBMLSpatial = any(strcmp(output_formats, 'sbml')) && (isempty(reaction_network_file) || any(strcmpi(reaction_network_file_ext, {'.xml'})));
        options.output.VCML.writeVCML = any(strcmp(output_formats, 'vcml')) && (isempty(reaction_network_file) || any(strcmpi(reaction_network_file_ext, {'', '.net', '.vcml'})));
        options.output.writeMCellMDL = any(strcmp(output_formats, 'mcell')) && (isempty(reaction_network_file) || any(strcmpi(reaction_network_file_ext, {'.net', '.mdl', '.mcell'})));

        options.VCML.downsampling = downsampling;

        prefix = 'img';
        if any(strcmpi(reaction_network_file_ext, {'.net', '.vcml'}))
            prefix = sprintf('%s_%s', prefix, reaction_network_file_name);
        elseif any(strcmpi(reaction_network_file_ext, {'.mdl', '.mcell'}))
            prefix = sprintf('%s_%s', prefix, strrep(reaction_network_file_name, '.*', ''));
        end
        
        if isempty(reaction_network_file)
        elseif strcmpi(reaction_network_file_ext, '.net')
            options.NET.filename = reaction_network_file;
            options.NET.units.length = 'm';
            options.NET.units.concentration = 'mol.m-3';
            options.NET.useImageAdjacency = false;
        elseif strcmpi(reaction_network_file_ext, '.vcml')
            options.VCML.input_filename = reaction_network_file;
            options.NET.useImageAdjacency = false;
        elseif any(strcmpi(reaction_network_file_ext, {'.mdl', '.mcell'}))
            options.MCellMDL.input_filename_pattern = reaction_network_file;
            options.NET.useImageAdjacency = false;
        else
            error('reaction network file of type ''%s'' not supported', reaction_network_file_ext);
        end
        
        
        
        prefix = sprintf('%s_%2.2f', prefix, downsampling);
        options.prefix = prefix;

        options.NET.translations = {'cell', 'CP'; 'nuc', 'NU'; 'nucleus', 'NU'; 'lamp2_mat_tfr_mat', 'EN'; 'CP_EC', 'PM'; 'CP_EN', 'EM'; 'CP_NU', 'NM'};
        options.VCML.translations = options.NET.translations;
        options.SBML.translations = options.NET.translations;

        options.VCML.default_time_step = 1e0;
        options.VCML.end_time = 4000;
        suffix = sprintf('%s_%04dsec', suffix, options.VCML.end_time);
        
        options.VCML.output_time_step = 100;
        options.VCML.max_time_step = 4;
        options.VCML.absolute_tolerance = 1e-8;
        options.VCML.relative_tolerance = 1e-8;

        options.MCellMDL.end_time = options.VCML.end_time;
        options.MCellMDL.default_time_step = options.VCML.default_time_step;
        options.MCellMDL.max_time_step = options.VCML.max_time_step;
        options.MCellMDL.output_time_step = options.VCML.output_time_step;

        options.MCellMDL.interaction_radius = 0.03;

        options.SBML.end_time = options.VCML.end_time;
        options.SBML.default_time_step = options.VCML.default_time_step;
        options.SBML.max_time_step = options.VCML.max_time_step;
        options.SBML.output_time_step = options.VCML.output_time_step;




        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        % Produce files for at least n_images_to_synthesize simulations
        seed_base = reaction_network_file_base_seed;
        seed_index = 0;
        mdl_pattern = [prefix, '*', suffix, filesep, 'cell1', filesep, 'cell.geometry.mdl'];
        mdl_count_function = @()length(dir(mdl_pattern));
        seeds_finished = containers.Map('KeyType', 'uint32', 'ValueType', 'logical');
        while mdl_count_function() < n_images_to_synthesize
            seed = seed_base + seed_index;
            seed_index = seed_index + 1;
            
            seeds_finished(seed) = false;
            [seed_prefix, seed_temporary_results, chunk_filename, can_start, final_name, final_exists, temp_name] = custom_chunk_start(seed, prefix, suffix);
            
            if ~can_start
                if final_exists
                    seeds_finished(seed) = true;
                end
                continue;
            end

            slml2img_start_time = tic;
            slml2img_start_cputime = cputime;
            
            options.seed = seed;
            options.prefix = seed_prefix;
            options.temporary_results = seed_temporary_results;
            
            options.numberOfGaussianObjects = n_gaussian_objects_function(seed);
            
            try
                state = rng( options.seed );
            catch err
                rand( 'seed', options.seed ); %#ok<RAND>
            end
            answer = slml2img( {combined_models}, options );

            slml2img_elapsed_time = toc(slml2img_start_time);
            slml2img_elapsed_cputime = cputime - slml2img_start_cputime;
            fprintf('\nslml2img for %s took %.3f s (%.3f s CPU time)\n\n', seed_prefix, slml2img_elapsed_time, slml2img_elapsed_cputime);
            
            if answer
                empty = [];
                save(final_name, 'empty');
                seeds_finished(seed) = true;
            end
            
            chunk_finish(chunk_filename);
        end


        if strcmpi(reaction_network_file_ext, '.vcml') && all(cell2mat(seeds_finished.values()))
            % Combine all VCML files into one for convenience
            
            combined_vcml_name = sprintf('%s_%010d_%04d%s', prefix, reaction_network_file_base_seed, n_images_to_synthesize, suffix)
            [combined_can_start, combined_final_name, combined_final_exists, combined_temp_name] = chunk_start(combined_vcml_name, '.vcml');
            
            if combined_can_start
                fprintf('\nCombining VCML files\n');
                combining_start_time = tic;
                combining_start_cputime = cputime;
                
                all_vcml = [];
                all_vcml_biomodel = [];
                for seed = seeds
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
    end
end



elapsed_time = toc(start_time);
elapsed_cputime = cputime - start_cputime;
fprintf('\n%s took %.3f s (%.3f s CPU time)\n\n', mfilename, elapsed_time, elapsed_cputime);

end%generate_simulation_instances_SarmaGhosh2012
