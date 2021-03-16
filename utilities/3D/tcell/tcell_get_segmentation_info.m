function [t_cell_info] = tcell_get_segmentation_info(t_cell_info)
  % Variables referenced more than once can be copied to local variables:
  master_script_options = t_cell_info.options;
  regions_location = t_cell_info.path_info.regions_location;
  segmentations_location = t_cell_info.path_info.segmentations_location;
  % condition_retrieval_function = t_cell_info.synapse_info.condition_retrieval_function;
  % sensors = t_cell_info.synapse_info.sensors;
  % number_sensors = t_cell_info.synapse_info.number_sensors;
  % sensor_retrieval_function = t_cell_info.synapse_info.sensor_retrieval_function;
  % number_runs = t_cell_info.synapse_info.number_runs;
  window_center_2d = t_cell_info.preprocessing_info.window_center_2d;
  window_size_2d = t_cell_info.preprocessing_info.window_size_2d;
  
  template_cell_radius = t_cell_info.template_info.template_cell_radius;
  background_subtraction_function = t_cell_info.preprocessing_info.background_subtraction_function;
  
  sensor_specific_segmentation_coarse_snake_options = master_script_options.sensor_specific_segmentation_coarse_snake_options;
  sensor_specific_segmentation_fine_snake_options = master_script_options.sensor_specific_segmentation_fine_snake_options;
  
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

  template_image = t_cell_info.template_info.template_image;

  window_size_2d = t_cell_info.preprocessing_info.window_size_2d;
  window_center_2d = t_cell_info.preprocessing_info.window_center_2d;
  cropped_size = t_cell_info.preprocessing_info.cropped_size;
    
  
  template_volume = sum(template_image(:));
  segmentation_volume_threshold = template_volume * master_script_options.segmentation_volume_minimum_ratio;
  t_cell_info.segmentation_info.template_volume = template_volume;
  t_cell_info.segmentation_info.segmentation_volume_threshold = segmentation_volume_threshold;
end