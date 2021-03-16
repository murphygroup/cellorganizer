function [segmentation_structure] = two_stage_snake_segmentation(preprocessing_structure, options)
  % Do snake segmentation twice; once coarsely and again to refine. Not very generic, made to be used in master_script_segment_windows.
  % 
  % 
  % 2014-04-25 tebuck: Copied from master_script_segment_windows.m.
  % 
  % 
  % 
  % 
  % 
  % % Dependencies:
  % % From the File Exchange: export_fig, pmkmp, toolbox_graph, toolbox_fast_marching, affine, Snake3D, Mesh_voxelisation, dirr, inhull
  % % From tebuck: modified version of Jieyue Li's HPA_lib, rasterize_mesh

  
  default_options = struct();
  % default_options.background_subtraction_function = @(x)x;
  default_options.background_subtraction_function = @(x)x - min(x(:));
  default_options.window_center_2d = [];
  default_options.seed_radius = [];
  default_options.cell_diameter = [];
  default_options.coarse_snake_options = struct;
  default_options.fine_snake_options = struct;
  
  
  if ~exist('options', 'var')
    options = default_options;
  else
    options_prior_to_processing = options;
    options = process_options_structure(default_options, options);
  end
  
  
  function [segmentation_structure] = segment_image(given_preprocessing_structure, given_options)
    synapse_locations_3d = given_preprocessing_structure.synapse_locations_3d;
    refined_synapse_locations_3d = given_preprocessing_structure.refined_synapse_locations_3d;
    image_to_segment = given_preprocessing_structure.image_to_segment;
    coherent_image_to_segment = given_preprocessing_structure.coherent_image_to_segment;
    coarse_white_line_image_to_segment = given_preprocessing_structure.coarse_white_line_image_to_segment ;
    fine_white_line_image_to_segment = given_preprocessing_structure.fine_white_line_image_to_segment;
    
    % Using parameters from/similar to r183 segment_test_cells_using_snakes.m and give the initial seed as a mesh:
    
    default_snake_options = struct(...
      'relative_convergence_criterion', 1e-2 ...
      , 'seed_subdivisions', 3 ...
      , 'Gamma', 0.5 ...
      , 'Iterations', 10 ...
      , 'Sigma1', 2 ...
      , 'Wline', -1e-2 ...
      , 'Wedge', 5e-2 ...
      , 'Sigma2', 2 ...
      , 'GIterations', ceil((options.cell_diameter / 2) * .4) ...
      , 'Sigma3', 3 ...
      , 'Alpha', 1.2 ...
      , 'Beta', 0.5 ...
      , 'Delta', 1e-2 ...
      , 'Kappa', 2e1 ...
      , 'Verbose', false ...
      , 'replicate_padding', true ...
      , 'save_last_frame', false ...
      , 'voxel_scales', [1, 1, 1/4] ...
      );
      % , 'GIterations', (ceil(template_cell_radius * 0.25) + 4) ...
    if ~exist('given_options', 'var')
      given_options = default_snake_options;
    else
      given_options = process_options_structure(default_snake_options, given_options);
    end
    
    segmentation_mesh = [];
    if isfield(given_options, 'segmentation_mesh')
      [segmentation_mesh] = snake_segmentation(double(fine_white_line_image_to_segment), given_options.segmentation_mesh, options.seed_radius, given_options);
    else
      [segmentation_mesh] = snake_segmentation(double(fine_white_line_image_to_segment), refined_synapse_locations_3d, options.seed_radius, given_options);
    end
    
    segmentation_structure = struct();
    segmentation_structure.segmentation_mesh = segmentation_mesh;
    segmentation_structure.given_options = given_options;
  end

  
  % Given a preprocessed image structure from two_stage_snake_segmentation_window_preprocessing, segment twice using snake_segmentation, once with seed points, another time with different options and starting with the output of the first segmentation:
  v = segment_image(preprocessing_structure, options.coarse_snake_options);
  v2 = segment_image(preprocessing_structure, setfield(options.fine_snake_options, 'segmentation_mesh', v.segmentation_mesh));
  segmentation_structure = v2.segmentation_mesh;
  
  
end

