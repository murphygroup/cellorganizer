
function [options] = tcell_process_options(options)
  % 2016-03-06 tebuck: Copied from master_script_process_options.m.
  % % Dependencies:
  % % From the File Exchange: export_fig, pmkmp, toolbox_graph, toolbox_fast_marching, affine, Snake3D, Mesh_voxelisation, dirr, inhull
  % % From tebuck: modified version of Jieyue Li's HPA_lib, rasterize_mesh
  % 02/23/2016 xruan: simplify code and delete unnecessary code
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Set defaults and process options structure:
  
  
  default_options = tcell_get_default_options();
  
  if ~exist('options', 'var')
    options = default_options;
  else
    options_prior_to_processing = options;
    options = process_options_structure(default_options, options);
  end

  
  % Allow shortcuts to frequently used data sets for convenience of development:
  
  % First create a base set of options for some all_conditions variants:
  % Use Excel files for all conditions:  
  % all_conditions_base_options.image_location = {'/projects/cellorganizer/xruan/tcell_project/images_organized/tcell_paper_images/'};
  % all_conditions_base_options.results_location = '/projects/cellorganizer/xruan/tcell_project/tcell_excel_results/';
  % all_conditions_base_options.use_segmentation_filtering = false;
  % all_conditions_base_options.use_segmentation_filtering = true;
  % all_conditions_base_options.skip_morphing_inside_template = true;
  % all_conditions_base_options.sensors_to_exclude = {};
  % all_conditions_base_options.sensors_to_include = {'Actin'};
  % all_conditions_base_options.timepoints_to_include = [0];
  % all_conditions_base_options.conditions_to_include = {'Full Stimulus'};
  % all_conditions_base_options.model_type_to_include = {'standardized_voxels'};
  
  % Ignore:
  % all_conditions_base_options.sensors_to_exclude{end + 1, 1} = 'Vav1';
  % Full Stimulus HS1 has poorly segmented high magnification runs that dominate the model:
  % xruan 02/23/2016 disable warning
  % warning('2014-11-19: Add rule for HS1 images! Or delete associated .csv files!')
  % xruan 08/04/2015 
  % include HS1
  % all_conditions_base_options.sensors_to_exclude{end + 1, 1} = 'HS1';
  % all_conditions_base_options.sensors_to_exclude{end + 1, 1} = 'Ezrin';
  % all_conditions_base_options.sensors_to_exclude{end + 1, 1} = 'Akt';
  % all_conditions_base_options.condition_name_abbreviations = {'B7 Blockade active Rac Cofilin - ', 'B7 ac Rac Cfl - '};
  % ll_conditions_base_options.condition_name_abbreviations = {'B7 Blockade active Rac Cofilin', 'B7BlkActiveRacCfl'; 'B7 Blockade', 'B7Blk'; 'Full Stimulus', 'FullStim'};
  
  
  % options.debug_use_debug_switch_values_in_master_script = true;
  if options.debug_use_debug_switch_values_in_master_script
    % For debugging, skip everything that is already (partially) done:
    options.skip_window_saving = true;
    options.skip_segmentation = true;
    options.skip_alignment = true;
    options.skip_morphing = true;
    % % Takes maybe a minute, run every so often:
    options.skip_morphing_inside_template = true;
    options.skip_model_building = true;
    options.skip_statistical_testing_across_conditions = true;
    % options.skip_hierarchical_clustering = true;
    options.skip_model_figure_generation = true;
    options.skip_model_video_saving = true;
    % options.skip_illustration_saving = true;
    % options.skip_create_blender_figures = true;
  end
  
  
  % options, pause
  
end

