function [t_cell_info, options] = tcell_build_models(t_cell_info, options)
% The function is called in img2model.m after parameterizing of the cells. 
% It uses t_cell_info and options as input and output. 
% 
% The basic idea is to take all standardized cells in, and build an average 
% model of the cells voxel by voxel. If there is input of different time points, 
% then it will build models for each time point. 
% 
% To specify the methods for the model type, it first calls tcell_get_model_type_info.m. 
% Then it builds the model following the algorithm of the model type. 
% 
% The model(s) is saved in the disk. Here the model is not the valid type of 
% model for CellOrganizer, and we need to further process the model.  
% 
% 2016-02-29: xruan: Copied from master_script_build_models.m.
% 
% 
% 
% 
% 
% % Dependencies:
% % From the File Exchange: export_fig, pmkmp, toolbox_graph, toolbox_fast_marching, affine, Snake3D, Mesh_voxelisation, dirr, inhull
% % From tebuck: modified version of Jieyue Li's HPA_lib, rasterize_mesh
%
% Author: Xiongtao Ruan
 

% Variables referenced more than once can be copied to local variables:

% get model type representation. 

t_cell_info = options.t_cell_info; 
[t_cell_info] = tcell_get_model_type_info(t_cell_info);

master_script_options = t_cell_info.options;
regions_location = t_cell_info.path_info.regions_location;
segmentations_location = t_cell_info.path_info.segmentations_location;
segmentations_filtered_location = t_cell_info.path_info.segmentations_filtered_location;
alignments_location = t_cell_info.path_info.alignments_location;
deformations_location = t_cell_info.path_info.deformations_location;
models_location = t_cell_info.path_info.models_location;

run_relative_path_list = t_cell_info.synapse_info.run_relative_path_list;
run_key_list = t_cell_info.synapse_info.run_key_list;
% run_data_list = t_cell_info.synapse_info.run_data_list;
run_images_list = t_cell_info.synapse_info.run_images_list;

timepoints_to_include = options.timepoints_to_include;
all_run_data = t_cell_info.synapse_info.all_run_data;
num_imgs = size(all_run_data, 1);
relative_time_run_list = t_cell_info.synapse_info.relative_time_run_list;
frame_tracking_list = t_cell_info.synapse_info.frame_tracking_list;
synapse_tracks_run_list = t_cell_info.synapse_info.synapse_tracks_run_list;

model_types = t_cell_info.model_type_info.model_types;
number_model_types = t_cell_info.model_type_info.number_model_types;
model_representation_functions = t_cell_info.model_type_info.model_representation_functions;
model_reconstruction_functions = t_cell_info.model_type_info.model_reconstruction_functions;
model_dimensionality_reduction_methods = t_cell_info.model_type_info.model_dimensionality_reduction_methods;
model_dimensionality_reduction_method_parameters = t_cell_info.model_type_info.model_dimensionality_reduction_method_parameters;
model_dimensionality_reduction_method_transformations = t_cell_info.model_type_info.model_dimensionality_reduction_method_transformations;

template_centroid = t_cell_info.template_info.template_centroid;
template_synapse = t_cell_info.template_info.template_synapse;
template_image = t_cell_info.template_info.template_image;

window_size_2d = t_cell_info.preprocessing_info.window_size_2d;
window_center_2d = t_cell_info.preprocessing_info.window_center_2d;
cropped_size = t_cell_info.preprocessing_info.cropped_size;
background_level_function = t_cell_info.preprocessing_info.background_level_function;

segmentation_rasterization_function = t_cell_info.segmentation_info.segmentation_rasterization_function;

template_volume = t_cell_info.segmentation_info.template_volume;
segmentation_volume_threshold = t_cell_info.segmentation_info.segmentation_volume_threshold;

% debug_print_model_filenames = false;
debug_print_model_filenames = true


% condition_sensor_combinations, number_condition_sensors, conditions, sensors, number_conditions, number_sensors
model_types, number_model_types
% number_all_relative_times

if master_script_options.skip_model_building

warning('>>>> HACK, options.skip_model_building true!')

else 
        
  % We want separate model files for each relative time to keep down memory usage and sooner have intermediate results:
  timepoints_to_include = t_cell_info.synapse_info.included_timepoints;
  for all_relative_time_index = 1 : numel(timepoints_to_include)
    relative_time = timepoints_to_include(all_relative_time_index);

    fprintf(' all_relative_time_index %d, relative_time % d\n', all_relative_time_index, relative_time)

    % We want separate model files for each model type (makes adding new model types easier):
    cur_run_data = all_run_data(relative_time_run_list == relative_time, :);
    cur_frame_index_list = frame_tracking_list(relative_time_run_list == relative_time, :);
    cur_synapse_run_list = synapse_tracks_run_list(relative_time_run_list == relative_time, :);

    for model_type_index = 1:number_model_types
      model_type = model_types{model_type_index};
      representation_function = model_representation_functions{model_type_index};
      reconstruction_function = model_reconstruction_functions{model_type_index};

      % This is the .mat file, not the figure:
      current_model_filename = sprintf('%sreltime_%d', master_script_options.model_prefix, relative_time);
      % current_model_filename = sprintf('%s_%s_frame%03d_reltime% 3d_%s', current_condition, current_sensor, all_relative_time_index, relative_time, model_type);
      % xruan 02/24/2016
      % [can_start, final_name, final_exists] = chunk_start(models_location, current_model_filename);
      [can_start, final_name, final_exists] = chunk_start_clean(models_location, current_model_filename);

      if debug_print_model_filenames, disp(sprintf('Model file ''%s'': can_start %d, final_exists %d\n', current_model_filename, can_start, final_exists)), end

      if ~can_start
        continue
      end

      % 02/25/2016 xruan add clean of old tmp files
      if ~can_start && ~final_exists
          chunck_lock_clean(models_location);
      end

      fprintf('tcell model_type_index %d, tcell model_type "%s"\n', model_type_index, model_type)

      current_count = 0;

      % Running/accumulated model:
      current_mean = [];
      current_standard_deviation = [];
      current_covariance = [];

      should_compute_covariance = false;

      % Batch model:
      % Preallocate as much space as might be necessary:
      current_number_cells = size(cur_run_data, 1);
      current_number_features = numel(representation_function(double(template_image)));
      current_data = zeros(current_number_cells, current_number_features);

      for i = 1 : current_number_cells  
          frame_index = cur_frame_index_list(i, :);
          current_synapse_annotations = cur_synapse_run_list(i, :);
          current_synapse_center_rounded = round([mean(current_synapse_annotations([1, 3])), mean(current_synapse_annotations([2, 4]))]);              
          cur_deformation_filename = [sprintf('run%d_cell%02d_frame%05d_synapse%05d,%05d', frame_index(2), frame_index(3), frame_index(1), current_synapse_center_rounded)];

          % Eventually, this should wait for all .mat files, and master_script_morph_segmentations.m should save .mat files with empty results if failed instead of leaving .tmp files. This way, running multiple jobs, others will quit when chunk_start here fails but the job that computes will get everything in the model files.

          if ~exist([deformations_location, cur_deformation_filename, '.mat'], 'file')
            if false && exist([segmentations_location, cur_deformation_filename, '.mat'], 'file')
              fprintf('   "%s" does not exist, but segmentation does!\n', [deformations_location, cur_deformation_filename, '.mat'])
              % beep, keyboard
            end
            if false && exist([segmentations_filtered_location, cell_deformation_filename, '.mat'], 'file')
              fprintf('   "%s" does not exist, but filtered segmentation does!\n', [deformations_location, cur_deformation_filename, '.mat'])
              % beep, keyboard
            end
            continue
          end

          % fprintf('   "%s" exists!\n', [deformations_location, cell_deformation_filename, '.mat'])

          try 
            a = load([deformations_location, cur_deformation_filename, '.mat']);
          catch 
              delete([deformations_location, cur_deformation_filename, '.mat'])
              disp(sprintf('unable to load %s', [deformations_location, cur_deformation_filename, '.mat']));
              return;
          end
          if isempty(fieldnames(a))
            continue
          end
          if ~isfield(a, 'cropped_raw_image_deformed')
              flag = 1;
          end
          current_window_image = a.current_window_image;     
          cropped_raw_image_deformed = a.cropped_raw_image_deformed;
          cropped_segmentation_image_deformed = a.cropped_segmentation_image_deformed;

          current_volume = sum(cropped_segmentation_image_deformed(:));
          is_segmentation_low_volume = current_volume <= segmentation_volume_threshold;

          % error('Correct if stmt')
          if is_segmentation_low_volume
            continue
          end
          if i > 50
              flag = 1;
          end
          % Subtract the original image's estimated background, then normalize total intensity to one:
          current_deformed_image = cropped_raw_image_deformed - background_level_function(current_window_image);
          current_deformed_image(current_deformed_image < 0) = 0;
          current_deformed_image(template_image < .5) = 0;
          current_deformed_image = current_deformed_image ./ sum(current_deformed_image(:));

          current_parameters = representation_function(current_deformed_image);

          % Batch model:
          current_count = current_count + 1;
          % current_data(current_count, :) = reshape(current_parameters, 1, []);
          current_data(current_count, :) = current_parameters;
          % keyboard              
      end

      % Batch model:
      current_data = current_data(1:current_count, :);
      current_mean = mean(current_data, 1);
      current_standard_deviation = std(current_data, [], 1);
      if should_compute_covariance
        current_covariance = cov(current_data);
      end

      % if any(isnan(current_mean))
        % beep, keyboard
      % end

      fprintf('   current_count %d\n', current_count)

      if false && current_count == 0
        % Debug info:
        beep, keyboard
      end

      if current_count > 0
        save(final_name, 'current_data', 'current_mean', 'current_standard_deviation', 'current_covariance', 'relative_time', 't_cell_info')
        % copy model to current directory
        % Leave the .tmp lock file so this work isn't repeated:
        chunk_finish(models_location, current_model_filename);
        system(sprintf('mv "%s" %s', final_name, options.results_location));
      end

      % Leave the .tmp lock file so this work isn't repeated:
      % chunk_finish(models_location, current_model_filename);

      current_count = 0;
      % pause
    end

  end

  % Compute dimensionality reduction transforms (this is required for later code, so maybe put "if master_script_options.skip_model_building" in this file, then use chunk_start for all of the above work?):
end

% beep, keyboard
% error('Implementation yet unfinished below this line!')
% warning('Implementation yet unfinished below this line!')

t_cell_info.model_info = struct();
t_cell_info.model_info.model_representation_functions = model_representation_functions;
t_cell_info.model_info.model_reconstruction_functions = model_reconstruction_functions;

options.t_cell_info = t_cell_info;
  
end  

