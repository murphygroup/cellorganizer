function [default_options] = tcell_get_default_options()
  % 2016-03-06 xruan: Copied from master_script_get_default_option.m.
  % 2013-07-18 tebuck: Adding what were previously variables in master_script.m.

  default_options.Wrapper_root = 'bhcho_bin_Tcell2/Wrapper/';
  default_options.image_location = '/images/Wuelfing_New_data/KTR CMU/';
  default_options.synapse_location = '/home/bhcho/bin/Tcell2/New_data_folders_only/KTR_CMU/';
  % default_options.synapse_location = 'bhcho_bin_Tcell2/New_data_folders_only/';
  default_options.results_location = './models/';
  default_options.tmp_results_location = './temp/';
  % default_options.annotation_location = '/projects/cellorganizer/xruan/tcell_project/tcell_annotation/';
  default_options.set_default_plot_options = true;
  default_options.use_profiling = false;
  default_options.named_options_set = 'none';
  default_options.model_prefix = '';
  
  % If this is true, in the future the alignment code should use both synapse points when available to precisely align cells:
  default_options.use_two_point_synapses = false;
  
  % List of sensors whose coordinate files should be completely ignored:
  default_options.run_keys_to_exclude = {};
  default_options.sensors_to_exclude = {};
  default_options.date_ranges_to_exclude = zeros(0, 6);
  
  % in default, there is only one channel in the 3D stacks
  default_options.image_channels = 1;
  
  % Filter cells that appear in fewer than this number of frames:
  default_options.minimum_number_frames_per_cell = 2;

  % Produce a text file listing the position of each cell in its first frame for manual evaluation:
  % default_options.create_cell_position_text_file = false;
  default_options.create_cell_position_text_file = true;

  % Sensor, really, not just proteins:
  % % This is a cell array where the first column is the name of a particular sensor and the second column one of the alternate names:
  % This is a cell array where each entry is a cell array of variant names of a particular sensor:
  default_options.protein_name_synonyms = {{'CPalpha1', 'CPa1'}};
  
  % N x 2 cell array, first column is full condition name, second column is abbreviated condition name:
  default_options.condition_name_abbreviations = {};
  
  default_options.possible_image_extensions = {'tif', 'stk'};
  
  default_options.seed_radius = 10;
  % This is approximate:
  default_options.cell_diameter = 23;
  % default_options.room_scale = 1;
  default_options.room_scale = 1.5;
  % default_options.room_scale = 2;
  
  default_options.segmentation_coarse_snake_options = struct(...
    'Gamma', 2e-1 ...
    , 'Sigma1', 1 ...
    , 'Wline', -5e0 ...
    , 'Wedge', 0 ...
    , 'Sigma2', 2 ...
    , 'Sigma3', 2 ...
    , 'Alpha', 0.15 ...
    , 'Beta', 0.1 ...
    , 'Delta', 0.1 ...
    , 'Kappa', 1 ...
    , 'Iterations', 240 ...
    , 'GIterations', 0 ...
    , 'Lambda', .95 ...
    );
    % 'seed_subdivisions', 3 ...
    % , 'replicate_padding', true ...
    % , 'Gamma', 2e-1 ...
  default_options.segmentation_fine_snake_options = struct(...
    'Gamma', 5 ...
    , 'Sigma1', 1 ...
    , 'Wline', -5 ...
    , 'Wedge', 0 ...
    , 'Sigma2', 1 ...
    , 'Sigma3', 1 ...
    , 'Alpha', 0.15 ...
    , 'Beta', 0.1 ...
    , 'Delta', 0.1 ...
    , 'Kappa', 10 ...
    , 'Iterations', 240 ...
    , 'GIterations', ceil((default_options.cell_diameter / 2) * .4) ...
    , 'Lambda', .75 ...
    );
    % , 'GIterations', ceil(template_cell_radius * .2 / .5) ...
  
  
  
  % Segmentation parameters specific to particular sensors:
  default_options.sensor_specific_segmentation_coarse_snake_options = {};
  default_options.sensor_specific_segmentation_fine_snake_options = {};
  % These were found automatically by optimization as done in master_script_optimize_segmentation_parameters.m:
  akt_segmentation_coarse_snake_options = struct(...
    'Gamma', 0.1544 ...
    , 'Sigma1', 1.1001 ...
    , 'Wline', -4.9818 ...
    , 'Wedge', -0.1292 ...
    , 'Sigma2', 2.0417 ...
    , 'Sigma3', 1.9867 ...
    , 'Alpha', 0.0731 ...
    , 'Beta', 0.0947 ...
    , 'Delta', 0.1010 ...
    , 'Kappa', 0.8835 ...
    , 'Iterations', 240.0001 ...
    , 'GIterations', 1.2546 ...
    , 'Lambda', 0.9837 ...
    );
  akt_segmentation_fine_snake_options = struct(...
    'Gamma', 5.1596 ...
    , 'Sigma1', 1.2808 ...
    , 'Wline', -5.1198 ...
    , 'Wedge', 0.0551 ...
    , 'Sigma2', 0.5026 ...
    , 'Sigma3', 0.8158 ...
    , 'Alpha', 0.2928 ...
    , 'Beta', 0.1673 ...
    , 'Delta', 0.0031 ...
    , 'Kappa', 3.7828 ...
    , 'Iterations', 240.0001 ...
    , 'GIterations', 5.0037 ...
    , 'Lambda', 0.7606 ...
    );
  default_options.sensor_specific_segmentation_coarse_snake_options(end + 1, 1:2) = {'Akt', akt_segmentation_coarse_snake_options};
  default_options.sensor_specific_segmentation_fine_snake_options(end + 1, 1:2) = {'Akt', akt_segmentation_fine_snake_options};
  
  
  default_options.segmentation_volume_minimum_ratio = 1e-1;
  
  
  default_options.use_segmentation_filtering = false;
  default_options.segmentation_filtering_interactive = false;
  % default_options.segmentation_filtering_interactive = true;
  
  default_options.align_method = 'original';
  default_options.landmark_temporal_smoothing_method = 'none';
  % default_options.landmark_temporal_smoothing_proportion = .25;
  default_options.landmark_temporal_smoothing_proportion = .5;
  default_options.landmark_temporal_smoothing_outlier_relative_error = .5;
  
  default_options.intensity_morphing_inside_template_method = 'none';
  
  default_options.confidence_level = .95;
  % default_options.nonparametric_manova_number_permutations = 1000;
  % default_options.nonparametric_manova_number_permutations = 10000;
  default_options.nonparametric_manova_number_permutations = 1e5;
  % default_options.rda_manova_number_permutations = 1000;
  % default_options.rda_manova_number_permutations = 10000;
  default_options.rda_manova_number_permutations = 1e5;
  
  % % default_options.videos_contrast_enhancement_method = '';
  % default_options.videos_contrast_enhancement_method = 'none';
  % default_options.videos_contrast_enhancement_method = 'log-scale';
  % default_options.videos_contrast_enhancement_method = 'exp-scale';
  default_options.videos_contrast_enhancement_method = 'threshold';
  default_options.videos_show_all_channels = false;
  % default_options.videos_show_all_channels = true;
  % default_options.videos_color_headings = false;
  default_options.videos_color_headings = true;
  % default_options.videos_scale = 1;
  default_options.videos_scale = 2;
  
  default_options.illustration_uses_histogram_equalization = false;
  % default_options.illustration_uses_histogram_equalization = true;
  default_options.illustration_uses_logarithmic_scale = false;
  % default_options.illustration_uses_logarithmic_scale = true;
  
  default_options.figure_dpi = 150;
  % default_options.figure_maximum_width = round(6.5 * default_options.figure_dpi);
  default_options.figure_maximum_width = round(6.5 * get(0, 'ScreenPixelsPerInch'));
  % % The figure using export_fig gets resized, so while it should be 975 pixels wide, it ends up at 860 (but this is figure-specific, use figure_height_scale instead):
  % default_options.figure_maximum_width = round(6.5 * get(0, 'ScreenPixelsPerInch') * 975 / 860);
  
  % Steps to skip:
  default_options.skip_window_saving = false;
  % default_options.skip_window_saving = true;
  default_options.skip_segmentation = false;
  % default_options.skip_segmentation = true;
  default_options.skip_segmentation_filtering = false;
  % default_options.skip_segmentation_filtering = true;
  default_options.skip_alignment = false;
  % default_options.skip_alignment = true;
  default_options.skip_morphing = false;
  % default_options.skip_morphing = true;
  default_options.skip_morphing_inside_template = false;
  % default_options.skip_morphing_inside_template = true;
  default_options.skip_model_building = false;
  % default_options.skip_model_building = true;
  % default_options.skip_per_run_model_building = false;
  % default_options.skip_per_run_model_building = true;
  default_options.skip_model_figure_generation = false;
  % default_options.skip_model_figure_generation = true;
  default_options.skip_model_video_saving = false;
  default_options.skip_model_video_saving = true;
  default_options.skip_illustration_saving = false;
  % default_options.skip_illustration_saving = true;
  default_options.skip_statistical_testing_across_conditions = false;
  default_options.skip_statistical_testing_across_conditions = true;
  default_options.skip_hierarchical_clustering = false;
  % default_options.skip_hierarchical_clustering = true;
  default_options.skip_create_blender_figures = false;
  default_options.skip_create_blender_figures = true;
  default_options.skip_postprocess = true;
  
  % options for synapse inference
  default_options.infer_synapses = false;
  default_options.infer_starting_time_point = false;
  
  % options for alignment adjustment
  default_options.adjust_one_point_alignment = false;
  default_options.one_point_alignment_adjust_model_filename = 'MVR_SinCos_ZX_around_synapse_cylinder_ZXY_250_092301_only_ZX.mat';
  default_options.infer_theta_y_method = 'method_5';
  default_options.infer_theta_z_separately = true;
  
  
  % Debugging switches:
  default_options.debug_synapse_file_processing_verbose = false;
  default_options.debug_model_reconstruction = false;
  % Produce windows for only one cell for each run (so maximum number of windows processed is (number_runs * number_all_relative_times):
  default_options.debug_produce_windows_one_cell_per_run = false;
  % The above option can also use, e.g., two cells per run:
  default_options.debug_produce_windows_number_cells_per_run = 1;
  % default_options.test_model_types = false;

  
  % Feb. 12, 2016
  default_options.use_WAVE2_manual_seg = false;
  
  default_options.debug_use_debug_switch_values_in_master_script = false;
  
  if false
    default_options.debug_produce_windows_one_cell_per_run = true;
  end
  
  if default_options.skip_postprocess
      default_options.skip_hierarchical_clustering = true;
  end
      
  
end
