function [t_cell_info] = tcell_segment_windows(t_cell_info, options) 
% The inputs are t_cell_info, options, the former one contains the general 
% information for the pipeline, and the latter contains the information for 
% the specific cell. The output is  t_cell_info. And the output result relative 
% to the cell is saved in the disk.
% 
% The idea is to use the cropped image in the disk, and perform a two-stage 
% snake segmentation. The first stage is the coarse stage which is intended 
% to find the approximate contour that fit the cell and the second stage is 
% the fine stage that uses the contour from the first stage to refine the 
% segmentation. The technique detail can be referred to Royal et al. And the 
% segmentation is saved as a polygon mesh in the disk. 
% 
% 2016-02-29 xruan: Copied from master_script_segment_windows_new.m.
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
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Segment windows corresponding to individual synapses:

  
  % Variables referenced more than once can be copied to local variables:
  master_script_options = t_cell_info.options;
  regions_location = t_cell_info.path_info.regions_location;
  segmentations_location = t_cell_info.path_info.segmentations_location;  
  
  window_center_2d = t_cell_info.preprocessing_info.window_center_2d;
  window_size_2d = t_cell_info.preprocessing_info.window_size_2d;

  image_name = options.image_name;
  synapse_annotation = options.synapse_annotation;
  relative_time = options.relative_time;
  frame_channel = options.frame_channel;
  frame_index = options.frame_index;
  run_index = frame_index(:, 2);  
  
  template_cell_radius = t_cell_info.template_info.template_cell_radius;
  background_subtraction_function = t_cell_info.preprocessing_info.background_subtraction_function;
  
  % Add things other code might need to segmentation_info:
  
  % rasterization_oversampling_scale = 1;
  rasterization_oversampling_scale = 2;
  % rasterization_oversampling_scale = 4;
  segmentation_rasterization_function = @(given_mesh, window_size)rasterize_mesh(given_mesh, struct('oversampling_scale', rasterization_oversampling_scale, 'cropping', [ones(3, 1), window_size']));
  % segmentation_rasterization_function = @(given_mesh)rasterize_mesh(given_mesh, struct('oversampling_scale', rasterization_oversampling_scale, 'cropping', [ones(3, 1), [window_size_2d, window_number_slices]']));
  t_cell_info.segmentation_info.rasterization_oversampling_scale = rasterization_oversampling_scale;
  t_cell_info.segmentation_info.segmentation_rasterization_function = segmentation_rasterization_function;
  t_cell_info.segmentation_info.window_preprocessing_function = @two_stage_snake_segmentation_window_preprocessing;
  t_cell_info.segmentation_info.two_stage_segmentation = @two_stage_snake_segmentation;
  
  
  % Segmentation-specific options:
  two_stage_segmentation_options = struct;
  two_stage_segmentation_options.window_center_2d = window_center_2d;
  two_stage_segmentation_options.seed_radius = master_script_options.seed_radius;
  two_stage_segmentation_options.cell_diameter = master_script_options.cell_diameter;
  two_stage_segmentation_options.background_subtraction_function = background_subtraction_function;
  two_stage_segmentation_options.coarse_snake_options = master_script_options.segmentation_coarse_snake_options;
  two_stage_segmentation_options.fine_snake_options = master_script_options.segmentation_fine_snake_options;
  
  
  if master_script_options.skip_segmentation
  
    warning('>>>> HACK, options.skip_segmentation true!')
  
  else
    
    % error('Update to scale voxels to cubical, use proper param adjustments to scale!')
    % warning('Update to scale voxels to cubical, use proper param adjustments to scale!')
    
    % fprintf_diary_cycle_function('Segmenting window images for run "%s"\n', run_key)
    printed_run_work_message = false;
    
    % Go through frames of each cell:
    % for current_cell_index = reshape(find(any(cell_index_frame_indices_uncomputed{run_index}, 2)), 1, [])
    current_synapse_annotations = options.synapse_annotation;
    if isempty(current_synapse_annotations)
        fprintf('%s annotation is empty\n', current_synapse_annotations);
        return;
    end
    current_synapse_center_rounded = round([mean(current_synapse_annotations([1, 3])), mean(current_synapse_annotations([2, 4]))]);
    current_relative_time = relative_time;
    cell_segmentation_filename = [sprintf('run%d_cell%02d_frame%05d_synapse%05d,%05d', frame_index(2), frame_index(3), frame_index(1), current_synapse_center_rounded)];
    cell_region_filename = [regions_location, cell_segmentation_filename];
    [can_start, final_name, final_exists] = chunk_start_clean(segmentations_location, cell_segmentation_filename);
    if final_exists
        fprintf('Segmentation %s already exists\n', cell_segmentation_filename);
        return;
    end
    
    while ~final_exists
        if ~can_start && ~final_exists
            out_status = chunk_lock_clean(segmentations_location);
            if out_status
                disp('Delete lock files whose jobs are not running!\n')
            end
        end

        try 
            a = load([cell_region_filename, '.mat']);
            current_window_image = a.current_window_image;            
        catch
            save(final_name, '')
            chunk_finish(segmentations_location, cell_segmentation_filename);
            break;
        end
        if ~printed_run_work_message
          fprintf('Segmenting window image for frame "%s"\n', cell_segmentation_filename)
          printed_run_work_message = true;
        end;
        % Segment the 3D fluorescent image:
        window_number_slices = size(current_window_image, 3);
        current_two_stage_segmentation_options = two_stage_segmentation_options;
        % current_sensor = sensor_retrieval_function({run_key});
        % relevant_coarse_snake_option_indices = strcmpi(sensor_specific_segmentation_coarse_snake_options(:, 1), current_sensor);
        % relevant_fine_snake_option_indices = strcmpi(sensor_specific_segmentation_fine_snake_options(:, 1), current_sensor);

        % xruan 07/16/2015
        % check if segment the right image
        current_two_stage_segmentation_options.current_synapse_annotations_transformed = a.current_synapse_annotations_transformed;
        % xruan 07/21/2015
        % judge one or two point annotation in two stage segmentaiton
        current_two_stage_segmentation_options.use_two_point_synapses = master_script_options.use_two_point_synapses;            

        % Assumes custom parameters will either exist or not exist together for coarse and fine:
%         switch sum(relevant_coarse_snake_option_indices)
%           case 0
%             ;
%           case 1
%             current_two_stage_segmentation_options.coarse_snake_options = sensor_specific_segmentation_coarse_snake_options{relevant_coarse_snake_option_indices, 2};
%             current_two_stage_segmentation_options.fine_snake_options = sensor_specific_segmentation_fine_snake_options{relevant_fine_snake_option_indices, 2};
%           otherwise
%             error('More than one relevant segmentation option set for sensor %s!', current_sensor)
%         end
        [segmentation_structure] = two_stage_snake_segmentation(two_stage_snake_segmentation_window_preprocessing(current_window_image, current_two_stage_segmentation_options), current_two_stage_segmentation_options);

        % window_filename = cell_region_filename;
        save(final_name, 'segmentation_structure', 'window_number_slices')
        chunk_finish(segmentations_location, cell_segmentation_filename);

        final_exists = true;     
    end
  end
end

