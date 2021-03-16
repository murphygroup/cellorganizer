function [] = segment_using_snakes(just_run_finished_compute_function)
  % 2012-11-14 tebuck: Copied from illustrate_bk_method_segmentation_process_in_test_cells.m.
  % 2013-02-21 tebuck: Adding image alignment along with segmentation mesh.
  % 2013-03-14 tebuck: Adding segmentation filter.
  % 2013-04-02 tebuck: Adding time series alignment using rasterize_mesh.
  % Dependencies:
  % From the File Exchange: export_fig, pmkmp, toolbox_graph, toolbox_fast_marching, affine, Snake3D, Mesh_voxelisation, dirr, inhull
  % From tebuck: modified version of Jieyue Li's HPA_lib, rasterize_mesh

  start_time = tic;

  filename = mfilename
  base_filename = [filename, '/']; 

  warning('off', 'MATLAB:MKDIR:DirectoryExists');
  warning('off', 'Images:initSize:adjustingMag');
  mkdir(filename)

  n = getenv('PBS_JOBID');
  n = regexprep(n, '[\r\n]', '');
  if isempty(n)
    n = '';
  end
  date_text = datestr(now(), 'yyyymmddHHMMSSFFF')
  diary([base_filename, 'log', date_text, ...
         '_', n, '.txt']);

  if ~exist('just_run_finished_compute_function', 'var')
    just_run_finished_compute_function = false;
  end

  %set(0, 'DefaultAxesFontName', 'Times')
  %set(0, 'DefaultLineLineWidth', 3, 'DefaultLineMarkerSize', 15)
  set(0, 'DefaultLineLineWidth', 2, 'DefaultLineMarkerSize', 12)
  set(0, 'DefaultAxesFontName', 'Helvetica')
  % set(0, 'DefaultAxesFontSize', 18)
  %set(0, 'DefaultAxesFontSize', 14)
  set(0, 'DefaultAxesFontSize', 12)
  % set(0, 'DefaultAxesFontSize', 10)
  %set(gcf, 'PaperPositionMode', 'auto');
  
  
  % Should use help to automatically derive further dependencies from this file's dependencies...
  addpath(genpath('pmkmp'))
  addpath(genpath('export_fig'))
  addpath(genpath('toolbox_graph'))
  addpath(genpath('toolbox_fast_marching'))
  % % addpath(genpath('toolbox_wavelet_meshes'))
  addpath(genpath('affine'))
  addpath(genpath('BasicSnake_version2f'))
  addpath(genpath('DataHash'))
  addpath(genpath('Mesh_voxelisation'))
  % addpath(genpath('GrTheory'))
  % addpath(genpath('convnfft'))
  addpath(genpath('dirr'))
  addpath(genpath('ordfilt3'))
  addpath(genpath('coherencefilter_version5b'))
  % addpath(genpath('TriangleRayIntersection'))
  addpath(genpath('inpolyhedron'))
  addpath(genpath('cellstructeq'))
  addpath(genpath('HPA_lib'))
  addpath(genpath('lowner'))
  % addpath(genpath('../4dcell/inhull'));
  addpath('../4dcell');


  use_profiling = false
  %use_profiling = true

  if use_profiling
    profile('-memory','on'); profile reset; profile on;
  end
  
  % Set random seed:
  s = RandStream('mt19937ar', 'Seed', 892734);
  RandStream.setDefaultStream(s);

  
  base_options.Wrapper_root = '/projects/cellorganizer/xruan/tcell_project/tcell/bhcho_bin_Tcell2/Wrapper/';
  base_options.image_location = '/images/Wuelfing_New_data/KTR CMU/';
  % Test cells:
  % base_options.synapse_location = '/home/bhcho/bin/Tcell2/New_data_folders_only_initial_set/';
  % All cells:
  base_options.synapse_location = '/home/bhcho/bin/Tcell2/New_data_folders_only/KTR_CMU/';
  % base_options.use_deconvolution = true;

  bk_segmentation_location = '/projects/cellorganizer/xruan/tcell_project/tcell/bhcho_bin_Tcell2/Results/Features_New_Initial_Set3_127/';
  bk_histeq_segmentation_location = '/projects/cellorganizer/xruan/tcell_project/tcell/bhcho_bin_Tcell2/Results/Features_New_Initial_Set3_127_histeq/';
  bk_imghisteq_segmentation_location = '/projects/cellorganizer/xruan/tcell_project/tcell/bhcho_bin_Tcell2/Results/Features_New_Initial_Set3_127_imghisteq/';
  chosen_segmentation_location = '/projects/cellorganizer/xruan/tcell_project/tcell/bhcho_bin_Tcell2/Results/Features_New_Initial_Set3_127_chose_regular_vs_imghisteq/';
  snake_segmentation_location = '/projects/cellorganizer/xruan/tcell_project/tcell/bhcho_bin_Tcell2/Results/Features_New_Initial_Set3_127_snake/';
  snake_filtered_segmentation_location = '/projects/cellorganizer/xruan/tcell_project/tcell/bhcho_bin_Tcell2/Results/Features_New_Initial_Set3_127_snake_filtered/';
  snake_region_location = '/projects/cellorganizer/xruan/tcell_project/tcell/bhcho_bin_Tcell2/Results/Features_New_Initial_Set3_127_snake_regions/';
  snake_filtered_region_location = '/projects/cellorganizer/xruan/tcell_project/tcell/bhcho_bin_Tcell2/Results/Features_New_Initial_Set3_127_snake_regions_filtered/';
  snake_aligned_segmentation_location = '/projects/cellorganizer/xruan/tcell_project/tcell/bhcho_bin_Tcell2/Results/Features_New_Initial_Set3_127_snake_aligned/';
  
  use_test_cells = false;
  fprintf('>>>> HACK\n'), use_test_cells = true
  if use_test_cells
    snake_segmentation_location = '/projects/cellorganizer/xruan/tcell_project/tcell/bhcho_bin_Tcell2/Results/Features_New_Initial_Set3_127_test_snake/';
    snake_filtered_segmentation_location = '/projects/cellorganizer/xruan/tcell_project/tcell/bhcho_bin_Tcell2/Results/Features_New_Initial_Set3_127_test_snake_filtered/';
    snake_aligned_segmentation_location = '/projects/cellorganizer/xruan/tcell_project/tcell/bhcho_bin_Tcell2/Results/Features_New_Initial_Set3_127_test_snake_aligned/';
    snake_region_location = '/projects/cellorganizer/xruan/tcell_project/tcell/bhcho_bin_Tcell2/Results/Features_New_Initial_Set3_127_test_snake_regions/';
    snake_filtered_region_location = '/projects/cellorganizer/xruan/tcell_project/tcell/bhcho_bin_Tcell2/Results/Features_New_Initial_Set3_127_test_snake_regions_filtered/';
    base_options.synapse_location = '/home/bhcho/bin/Tcell2/New_data_folders_only_initial_set/';;
  end
  
  seed_radius = 10;
  % seed_radius = 15;
  % seed_radius = 20;
  % seed_radius = 25;
  % seed_radius = 30;
  
  % This is approximate:
  cell_radius = 23;
  
  % room = seed_radius * 2;
  % room = seed_radius * 4;
  room = cell_radius * 2;
  % room = cell_radius * 3;
  
  % fmin_use_single_cell = false
  fmin_use_single_cell = true
  fmin_use_single_image = false
  % fmin_use_single_image = true
  optimization_running = false
  
  function [preprocessing_structure] = compute_preprocessed_images(compute_function_arguments)
    if size(compute_function_arguments.frame_synapse_locations, 1) == 0
      preprocessing_structure = struct();
      return
    end
    

    % whos

    synapse_locations_3d = compute_function_arguments.frame_synapse_locations;
    % Z coordinate is not given (it is T/the frame number relative to synapse formation):
    synapse_locations_3d(:, 3) = round(size(compute_function_arguments.frame_image, 3) / 2);
    
    % Remove synapses outside the image:
    bad_value_locations = false(size(synapse_locations_3d));
    bad_value_locations(synapse_locations_3d < 1) = true;
    bad_value_locations(synapse_locations_3d(:, 1) > size(compute_function_arguments.frame_image, 2), 1) = true;
    bad_value_locations(synapse_locations_3d(:, 2) > size(compute_function_arguments.frame_image, 1), 2) = true;
    bad_value_locations(synapse_locations_3d(:, 3) > size(compute_function_arguments.frame_image, 3), 3) = true;
    bad_synapses = max(bad_value_locations, [], 2);
    synapse_locations_3d = synapse_locations_3d(~bad_synapses, :);

    if fmin_use_single_cell && optimization_running
      if size(synapse_locations_3d, 1) > 1
        synapse_locations_3d = synapse_locations_3d(1, :);
      end
    end
    
    mean_synapse_locations = mean(synapse_locations_3d(:, 1:2), 1);
    number_synapses = size(synapse_locations_3d, 1);
    
    image_to_segment = compute_function_arguments.frame_image;
    % Background subtraction (prevents edges at image boundaries):
    background_level = mode(image_to_segment(image_to_segment < mean(image_to_segment(:))));
    image_to_segment = image_to_segment - background_level;
    image_to_segment(image_to_segment < 0) = 0;
    % Scaling for consistency in force computations:
    image_to_segment = image_to_segment ./ prctile(image_to_segment(:), 99.9);
    % Hopefully makes gradients at higher intensities less influential:
    % image_to_segment = histeq(image_to_segment, 1024);
    image_to_segment = reshape(histeq(reshape_2d(image_to_segment), 64), size(image_to_segment));
    
    
    dic_image = compute_function_arguments.frame_dic_image;
    

    % Improve seed locations using nearby maxima:
    % Simple local gradient descent:
    % sigma = cell_radius * .5;
    % Prevent moving to adjacent cells by making it smaller:
    sigma = cell_radius * .2;
    % sigma = cell_radius * .5;
    % % Filter bad synapses that are in background and far from cells even after refinement:
    % high_intensity_synapses = interp3(smoothed_image_to_segment, synapse_locations_3d(:, 1), synapse_locations_3d(:, 2), synapse_locations_3d(:, 3)) > .1
    % synapse_locations_3d = synapse_locations_3d(high_intensity_synapses, :);
    % Skip this image if there are none:
    % clear smoothed_image_to_segment gx gy gz
    if size(synapse_locations_3d, 1) == 0
      preprocessing_structure = struct();
      return
    end

    
    preprocessing_structure.number_synapses = number_synapses;
    preprocessing_structure.synapse_locations_3d = cell(number_synapses, 1);
    preprocessing_structure.refined_synapse_locations_3d = cell(number_synapses, 1);
    preprocessing_structure.image_to_segment = cell(number_synapses, 1);
    preprocessing_structure.raw_image = cell(number_synapses, 1);
    preprocessing_structure.dic_image = cell(number_synapses, 1);
    preprocessing_structure.x_range = cell(number_synapses, 1);
    preprocessing_structure.y_range = cell(number_synapses, 1);
    preprocessing_structure.z_range = cell(number_synapses, 1);
    preprocessing_structure.coherent_image_to_segment2 = cell(number_synapses, 1);
    preprocessing_structure.fine_white_line_image_to_segment = cell(number_synapses, 1);
    preprocessing_structure.coarse_white_line_image_to_segment = cell(number_synapses, 1);
    for cell_index = 1:number_synapses
      % this is plot x, same as original synapse coordinate x (i.e., horizontal):
      x_range = max(min(floor(synapse_locations_3d(cell_index, 1) - room), size(image_to_segment, 2)), 1):max(min(ceil(synapse_locations_3d(cell_index, 1) + room), size(image_to_segment, 2)), 1);
      y_range = max(min(floor(synapse_locations_3d(cell_index, 2) - room), size(image_to_segment, 1)), 1):max(min(ceil(synapse_locations_3d(cell_index, 2) + room), size(image_to_segment, 1)), 1);
      % z_range = max(min(floor(synapse_locations_3d(cell_index, 3) - room), size(image_to_segment, 3)), 1):max(min(ceil(synapse_locations_3d(cell_index, 3) + room), size(image_to_segment, 3)), 1);
      z_range = 1:size(image_to_segment, 3);
      % z_range

      preprocessing_structure.x_range{cell_index} = x_range;
      preprocessing_structure.y_range{cell_index} = y_range;
      preprocessing_structure.z_range{cell_index} = z_range;
      preprocessing_structure.image_to_segment{cell_index} = image_to_segment(y_range, x_range, z_range);
      preprocessing_structure.raw_image{cell_index} = compute_function_arguments.frame_image(y_range, x_range, z_range);
      preprocessing_structure.dic_image{cell_index} = dic_image(y_range, x_range);
      
      % Faster to do this individually:
      preprocessing_structure.coherent_image_to_segment2{cell_index} = CoherenceFilter(preprocessing_structure.image_to_segment{cell_index}, struct('T', 50, 'dt', 2, 'Scheme', 'R', 'verbose', 'none'));
      [gx, gy, gz] = smooth_gradient(preprocessing_structure.coherent_image_to_segment2{cell_index}, 'scharr5', false, true);
      preprocessing_structure.fine_white_line_image_to_segment{cell_index} = sqrt(gx.^2 + gy.^2 + gz.^2);
      preprocessing_structure.coarse_white_line_image_to_segment{cell_index} = smooth3(preprocessing_structure.fine_white_line_image_to_segment{cell_index}, 'gaussian', 7, sqrt(3)).^2;
      
      preprocessing_structure.synapse_locations_3d{cell_index} = synapse_locations_3d(cell_index, :) - [x_range(1), y_range(1), z_range(1)] + 1;
      smoothed_image_to_segment = convnfft_fast(preprocessing_structure.image_to_segment{cell_index}, gaussian_kernel(sigma));
      preprocessing_structure.refined_synapse_locations_3d{cell_index} = climb_image(smoothed_image_to_segment, preprocessing_structure.synapse_locations_3d{cell_index});
    end
    
    % imshow(reshape_contrast(preprocessing_structure.coarse_white_line_image_to_segment, 3)), pause
    % And what about using the second derivatives to make nearby/touching surfaces unfavorable?
    
    % whos preprocessing_structure
  end

  
  function [segmentation_structure] = two_stage_segmentation(argument_structure, snake_options)
    % fprintf('>>>>>>>> HACK\n'), snake_options{1}.Verbose = true; snake_options{2}.Verbose = true;
    a = argument_structure;
    v = segment_image(a, snake_options{1});
    % fprintf('>>>>>>>> HACK, turn second stage back on!!!!!\n')
    if length(fieldnames(v)) == 0
      % keyboard
      v.argument_structure = rmfield(v.argument_structure, 'frame_image');
      v.argument_structure = rmfield(v.argument_structure, 'frame_dic_image');
      segmentation_structure = v;
      return
    end
    second_snake_options = snake_options{2};
    a2 = a;
    a2.previous_segmentation_mesh = v.segmentation_mesh;
    v = segment_image(a2, snake_options{2});
    % v = rmfield(v, 'argument_structure');
    % v.argument_structure
    % keyboard
    v.argument_structure = rmfield(v.argument_structure, 'frame_image');
    v.argument_structure = rmfield(v.argument_structure, 'frame_dic_image');
    segmentation_structure = v;
  end

  function [segmentation_structure] = segment_image(argument_structure, snake_options)
    % argument_structure
    if size(argument_structure.frame_synapse_locations, 1) == 0
      segmentation_structure = struct();
      return
    end
    % file_id = argument_structure.file_id;
    % fprintf([file_id.filename, '\n'])

    [preprocessing_structure] = compute_preprocessed_images(argument_structure);
    if length(fieldnames(preprocessing_structure)) == 0
      segmentation_structure = struct();
      return
    end

    number_synapses = preprocessing_structure.number_synapses;
    synapse_locations_3d = preprocessing_structure.synapse_locations_3d;
    refined_synapse_locations_3d = preprocessing_structure.refined_synapse_locations_3d;
    dic_image = preprocessing_structure.dic_image;
    x_range = preprocessing_structure.x_range;
    y_range = preprocessing_structure.y_range;
    z_range = preprocessing_structure.z_range;
    image_to_segment = preprocessing_structure.image_to_segment;
    coherent_image_to_segment2 = preprocessing_structure.coherent_image_to_segment2;
    coarse_white_line_image_to_segment = preprocessing_structure.coarse_white_line_image_to_segment ;
    fine_white_line_image_to_segment = preprocessing_structure.fine_white_line_image_to_segment;
    
    
    % segmentation_structure = segmentation_structure;

    
    % Using parameters from/similar to r183 segment_test_cells_using_snakes.m and give the initial seed as a mesh:

    % This is approximate:
    default_template_options = get_default_template_options();
    % cell_radius = 23;
    cell_radius = (default_template_options.yr + default_template_options.zr) / 2;

    default_snake_options = struct(...
      'relative_convergence_criterion', 1e-2 ...
      , 'seed_subdivisions', 3 ...
      , 'Gamma', 2e-1 * 2.5 ...
      , 'Iterations', 10 * 1 ...
      , 'Sigma1', 2 ...
      , 'Wline', -1e-2 * 1 ...
      , 'Wedge', 5e-2 * 1 ...
      , 'Sigma2', 2 ...
      , 'GIterations', (ceil(cell_radius * 0.25) + 4) ...
      , 'Sigma3', 3 * 1 ...
      , 'Alpha', 6e-1 * 2 + 0.2 * 0 ...
      , 'Beta', 1e0 * .5 + 0.2 * 0 ...
      , 'Delta', 1e-2 * 1 + 0.1 * 0 ...
      , 'Kappa', 2e1 * 1e0 + 2 * 0 ...
      , 'Verbose', false...
      , 'replicate_padding', true ...
      , 'save_last_frame', false ...
      , 'voxel_scales', [1, 1, 1/4] ...
      );
    if ~exist('snake_options', 'var')
      snake_options = default_snake_options;
    else
      snake_options = process_options_structure(default_snake_options, snake_options);
    end
    
    segmentation_mesh = cell(number_synapses, 1);
    for cell_index = 1:number_synapses
      if isfield(argument_structure, 'previous_segmentation_mesh')
        % [segmentation_mesh] = snake_segmentation(double(fine_black_line_image_to_segment), argument_structure.previous_segmentation_mesh, [], snake_options);
        [segmentation_mesh{cell_index}] = snake_segmentation(double(fine_white_line_image_to_segment{cell_index}), argument_structure.previous_segmentation_mesh{cell_index}, seed_radius, snake_options);
      else
        % [segmentation_mesh] = snake_segmentation(double(coarse_black_line_image_to_segment), refined_synapse_locations_3d, seed_radius, snake_options);
        % [segmentation_mesh] = snake_segmentation(double(coarse_white_line_image_to_segment), refined_synapse_locations_3d, seed_radius, snake_options);
        [segmentation_mesh{cell_index}] = snake_segmentation(double(fine_white_line_image_to_segment{cell_index}), refined_synapse_locations_3d{cell_index}, seed_radius, snake_options);
      end
    end
    
    segmentation_structure = struct();
    segmentation_structure.argument_structure = argument_structure;
    segmentation_structure.segmentation_mesh = segmentation_mesh;
    segmentation_structure.snake_options = snake_options;
    % Need this information for further mesh processing, e.g., synapse identification and alignment:
    segmentation_structure.number_synapses = number_synapses;
    segmentation_structure.synapse_locations_3d = synapse_locations_3d;
    segmentation_structure.refined_synapse_locations_3d = refined_synapse_locations_3d;
    segmentation_structure.dic_image = dic_image;
    segmentation_structure.x_range = x_range;
    segmentation_structure.y_range = y_range;
    segmentation_structure.z_range = z_range;

  end


  
  % function [] = segmentation_finished_compute_function(compute_function_arguments, is_interactive)
  function [plot_success] = segmentation_finished_compute_function(argument_structure, is_interactive, subplot_indices)
    % fprintf('In segmentation_finished_compute_function!\n')
    plot_success = false;
    
    % bk_filename = strrep(argument_structure.relative_path, '/', '_')
    fixed_relative_path = regexprep(argument_structure.relative_path, '/+', '/');
    fixed_relative_path = fixed_relative_path(find(fixed_relative_path == '/', 1, 'first') + 1:find(fixed_relative_path == '/', 1, 'last') - 1);
    % keyboard
    % pause

    % argument_structure
    if size(argument_structure.frame_synapse_locations, 1) == 0
      plot_success = false;
      return
    end
    
    % Check if any test cells are in this image:
    has_bk_segmentations = false;
    for cell_index = 1:size(argument_structure.frame_synapse_locations, 1)
      bk_filename = [strrep(fixed_relative_path, '/', '_'), num2str(argument_structure.frame_synapse_cell_ids(cell_index), '_%d'), num2str(argument_structure.frame_synapse_locations(cell_index, 3), '_%d'), num2str(argument_structure.frame_index, '_%d')];
      individual_bk_segmentation_path = [bk_segmentation_location, filesep, bk_filename, '.mat'];
      if exist(individual_bk_segmentation_path, 'file')
        has_bk_segmentations = true
        break
      end
    end
    if ~has_bk_segmentations
      plot_success = false;
      return
    end
    

    [preprocessing_structure] = compute_preprocessed_images(argument_structure);
    if length(fieldnames(preprocessing_structure)) == 0
      plot_success = false;
      return
    end

    % Plot the mesh segmentation:
    colors = pmkmp(256, 'CubicL');

    number_synapses = preprocessing_structure.number_synapses;
    % segmentation_mesh = cell(number_synapses, 1);
    for cell_index = 1:number_synapses
      bk_filename = [strrep(fixed_relative_path, '/', '_'), num2str(argument_structure.frame_synapse_cell_ids(cell_index), '_%d'), num2str(argument_structure.frame_synapse_locations(cell_index, 3), '_%d'), num2str(argument_structure.frame_index, '_%d')];
      individual_bk_segmentation_path = [bk_segmentation_location, filesep, bk_filename, '.mat'];
      if ~exist(individual_bk_segmentation_path, 'file')
        continue
      else
        fprintf(['@@@@ ', individual_bk_segmentation_path, ' exists!\n'])
        a = load(individual_bk_segmentation_path);
        % a.Best_extra_data
        if ~isfield(a.Best_extra_data, 'final_segmentation_without_z_centering')
          continue
        end
        bk_isosurface = isosurface(a.Best_extra_data.final_segmentation_without_z_centering, .5);
        bk_angle = a.Best_extra_data.degree_medial_axis + a.Best_extra_data.degree_synapse_plane_final;
      end
      % keyboard

      synapse_locations_3d = preprocessing_structure.synapse_locations_3d{cell_index};
      refined_synapse_locations_3d = preprocessing_structure.refined_synapse_locations_3d{cell_index};
      dic_image = preprocessing_structure.dic_image{cell_index};
      x_range = preprocessing_structure.x_range{cell_index};
      y_range = preprocessing_structure.y_range{cell_index};
      z_range = preprocessing_structure.z_range{cell_index};
      image_to_segment = preprocessing_structure.image_to_segment{cell_index};
      coherent_image_to_segment2 = preprocessing_structure.coherent_image_to_segment2{cell_index};
      coarse_white_line_image_to_segment = preprocessing_structure.coarse_white_line_image_to_segment{cell_index};
      fine_white_line_image_to_segment = preprocessing_structure.fine_white_line_image_to_segment{cell_index};

      chosen_image_to_view = coherent_image_to_segment2;
      % chosen_image_to_view = hessian_image_to_segment;
      % chosen_image_to_view = hessian_image_to_segment;

      % fprintf('here!\n')
      % keyboard
      segmentation_mesh = argument_structure.result.segmentation_mesh{cell_index};
      

      if ~exist('is_interactive', 'var')
        is_interactive = false;
        subplot_indices = [1, 1];
      end
        

 
      % Prepare images to be plotted to right of mesh segmentation plot:
      segmentation_result = segmentation_mesh;
      segmentation_result.vertices(:, 1) = segmentation_result.vertices(:, 1) + x_range(1) - 1;
      segmentation_result.vertices(:, 2) = segmentation_result.vertices(:, 2) + y_range(1) - 1;
      % segmentation_image = voxelise(x_range, y_range, z_range, segmentation_mesh);
      segmentation_image = voxelise(x_range, y_range, z_range, segmentation_result);
      segmentation_image = permute(segmentation_image, [2, 1, 3]);
      % images_to_plot_with_synapses = {};

      % images_to_plot_with_synapses{end+1} = mean((chosen_image_to_view .* ~segmentation_image) + (min(chosen_image_to_view(:)) .* segmentation_image), 3);

      number_rows = 2;
      number_columns = 4;
      % number_columns = 5;
      number_additional_columns = 1;
      % number_additional_columns = 2;
      if is_interactive
        % figure
      else
        clf
      end
      % scale = 2;
      scale = 3;
      set(gcf, 'Position', [100, 100, 100 * number_columns * scale, 100 * number_rows * scale])
      % set(gcf,'render','opengl')
      % set(gcf, 'Position', [100, 100, 400 * number_columns * scale, 100 * number_rows * scale])
      if ~is_interactive
        % fprintf('>>>>>>>> HACK\n')
        set(gcf, 'Visible', 'off');
      end
    
    

      % Plot images to right of mesh segmentation plot with synapses marked:
      % for subplot_index = 1:length(images_to_plot_with_synapses)
      % Last ones have the meshes overlaid:
    
      % subplot(number_rows, number_columns, 1:number_columns - 1)
      subplot(number_rows, number_columns, [(1:number_columns - number_additional_columns), (1:number_columns - number_additional_columns) + number_columns])
      chosen_image_to_view_with_synapse = chosen_image_to_view;
      % Highlight the synapses:
      rounded_synapse_locations_3d = round(synapse_locations_3d);
      % synapse_radius = 1;
      synapse_radius = 2;
      for synapse_index = 1:size(rounded_synapse_locations_3d, 1)
        chosen_image_to_view_with_synapse(rounded_synapse_locations_3d(synapse_index, 2) + (-synapse_radius:synapse_radius), rounded_synapse_locations_3d(synapse_index, 1) + (-synapse_radius:synapse_radius), rounded_synapse_locations_3d(synapse_index, 3)) = min(chosen_image_to_view(:));
        chosen_image_to_view_with_synapse(rounded_synapse_locations_3d(synapse_index, 2), rounded_synapse_locations_3d(synapse_index, 1), rounded_synapse_locations_3d(synapse_index, 3)) = max(chosen_image_to_view(:));
      end
      % whos chosen_image_to_view chosen_image_to_view_with_synapse
      
      minimum_level = min(chosen_image_to_view(:));
      background_level = prctile(chosen_image_to_view(:), 1);
      
      % Rotate by BK's segmentation angle:
      segmentation_image = imrotate(segmentation_image, bk_angle, 'bilinear', 'crop') >= .5;
      chosen_image_to_view = imrotate(chosen_image_to_view, bk_angle, 'bilinear', 'crop');
      chosen_image_to_view_with_synapse = imrotate(chosen_image_to_view_with_synapse, bk_angle, 'bilinear', 'crop');
      segmentation_mesh.vertices = [segmentation_mesh.vertices, ones(size(segmentation_mesh.vertices, 1), 1)];
      segmentation_mesh.vertices = segmentation_mesh.vertices * rotation_matrix3h(3, bk_angle);
      segmentation_mesh.vertices = segmentation_mesh.vertices(:, 1:3);
      
      chosen_image_to_view_with_synapse(chosen_image_to_view_with_synapse < minimum_level) = minimum_level;
      chosen_image_to_view(chosen_image_to_view < minimum_level) = minimum_level;
      
      number_reshape_rows = 4; imagesc(reshape_contrast([chosen_image_to_view_with_synapse; (chosen_image_to_view .* ~segmentation_image) + (background_level .* segmentation_image)], number_reshape_rows))
      % colormap(colors)
      colormap(gray)
      axis equal
      axis tight
      title({'Row 1: Coherent image', 'Row 2: Same without segmentation'})
      
      % subplot(number_rows, number_columns, number_rows * number_columns)
      % subplot(number_rows, number_columns, number_columns)
      possible_rows = 1:number_rows;
      % possible_rows = 1;
      possible_columns = number_columns - number_additional_columns + 1:number_columns;
      [temp_columns, temp_rows] = meshgrid(possible_columns, possible_rows);
      final_subplot_indices = temp_columns + number_columns .* (temp_rows - 1);
      final_subplot_indices = final_subplot_indices(:);
      subplot(number_rows, number_columns, final_subplot_indices)
      if ~is_interactive
        segmentation_mesh2 = segmentation_mesh;
        segmentation_mesh2.vertices = segmentation_mesh2.vertices - repmat(min(segmentation_mesh2.vertices, [], 1), size(segmentation_mesh2.vertices, 1), 1);
        h1 = patch(segmentation_mesh2, 'FaceColor', [1, .4, .1], 'EdgeColor', 'none');
        hold on
        bk_isosurface2 = bk_isosurface;
        bk_isosurface2.vertices = bk_isosurface2.vertices + repmat(-min(bk_isosurface2.vertices, [], 1) + max(segmentation_mesh2.vertices, [], 1) .* [1, 0, 0], size(bk_isosurface2.vertices, 1), 1);
        h2 = patch(bk_isosurface2, 'FaceColor', [.1, .4, 1], 'EdgeColor', 'none');
        hold off
        camlight, lighting phong
        axis equal
        axis tight
        legend([h1, h2], {'Snakes (r262)', 'BK'}, 'Location', 'SouthOutside')
        set(gca, 'ydir', 'reverse')
      end
      % subplot(number_rows, number_columns, [(number_rows - 1) * (number_columns - number_additional_columns + 1:number_columns), number_rows * number_columns])


      drawnow;

      if ~is_interactive
        % image_filename = [base_filename, filesep, bk_filename, '_original_refined'];
        image_filename = [base_filename, filesep, bk_filename, '_original_refined', num2str(cell_index, '_synapse%03d')];
        ppi = 90;
        % ppi = 150;
        % ppi = 300;
        export_fig(image_filename, '-opengl', '-png', '-a1', '-nocrop', ['-r', num2str(ppi)])
        % export_fig(image_filename, '-zbuffer', '-png', '-a1', '-nocrop')
        % export_fig(image_filename, '-painters', '-png', '-a1', '-nocrop')
      end
    end
    
    plot_success = true;
  end
  
  
  
  
  
  
  default_template_options = get_default_template_options();
  cell_radius = (default_template_options.yr + default_template_options.zr) / 2;
  % Now try to find parameters for File Exchange version of SnakeMoveIteration3D.m:
  initial_snake_options = struct(...
    'Gamma', 1 * 1e-1 ...
    , 'Sigma1', 2 ...
    , 'Wline', 1e-2 * 0 ...
    , 'Wedge', 5e-2 * 5 ...
    , 'Sigma2', 2 ...
    , 'Sigma3', 3 ...
    , 'Alpha', 6e-1 * .5 ...
    , 'Beta', 1 * .5 ...
    , 'Delta', 1e-2  * .5 ...
    , 'Kappa', 2e1 ...
    , 'Iterations', 15 ...
    , 'GIterations', (ceil(cell_radius) + 4) * 0 ...
    , 'Lambda', .8 ...
    );
  % integer_options = {'Iterations', 'GIterations'};
  snake_options_field_names = fieldnames(initial_snake_options);
  % snake_options_limits = struct(...
    % 'Gamma', [.5, 2] ...
    % , 'Sigma1', [.5, 5] ...
    % , 'Wline', [-1e1, 1e1] ...
    % , 'Wedge', [0, 1e1] ...
    % , 'Sigma2', [.5, 5] ...
    % , 'Sigma3', [.5, 5] ...
    % , 'Alpha', [1e-1, 1e1] ...
    % , 'Beta', [1e-1, 1e1] ...
    % , 'Delta', [1e-1, 1e1] ...
    % , 'Kappa', [1e-1, 1e1] ...
    % , 'Iterations', [15, 15] ...
    % , 'GIterations', [0, ceil(cell_radius) * 2] ...
    % , 'Lambda', [.01, .99] ...
    % );
  % snake_options_limits = struct(...
    % 'Gamma', [1e-2, 2e-1] ...
    % , 'Sigma1', [.5, 5] ...
    % , 'Wline', [-1e1, 1e1] ...
    % , 'Wedge', [-1e1, 1e1] ...
    % , 'Sigma2', [.5, 5] ...
    % , 'Sigma3', [.5, 5] ...
    % , 'Alpha', [1e-2, 1e1] ...
    % , 'Beta', [1e-2, 1e1] ...
    % , 'Delta', [1e-3, 1e1] ...
    % , 'Kappa', [1e-1, 1e2] ...
    % , 'Iterations', [15, 15] ...
    % , 'GIterations', [0, cell_radius * .5] ...
    % , 'Lambda', [.01, .99] ...
    % );
  % To get GIterations involved:
  snake_options_limits = struct(...
    'Gamma', [.1, 10] ...
    , 'Sigma1', [.5, 5] ...
    , 'Wline', [-20, 0] ...
    , 'Wedge', [-5, 5] ...
    , 'Sigma2', [.5, 5] ...
    , 'Sigma3', [.5, 5] ...
    , 'Alpha', [0, 5] ...
    , 'Beta', [0, 5] ...
    , 'Delta', [0, 5] ...
    , 'Kappa', [0.5, 5] ...
    , 'Iterations', [100, 500] ...
    , 'GIterations', [0, ceil(cell_radius * .5)] ...
    , 'Lambda', [.5, .99] ...
    );
  snake_options_to_use = cell2mat(struct2cell(snake_options_limits));
  snake_options_to_use = snake_options_to_use(:, 1) < snake_options_to_use(:, 2);
  snake_options_integer_valued = struct(...
    'Gamma', false ...
    , 'Sigma1', false ...
    , 'Wline', false ...
    , 'Wedge', false ...
    , 'Sigma2', false ...
    , 'Sigma3', false ...
    , 'Alpha', false ...
    , 'Beta', false ...
    , 'Delta', false ...
    , 'Kappa', false ...
    , 'Iterations', true ...
    , 'GIterations', true ...
    , 'Lambda', false ...
    );
  snake_options_integer_valued = cell2mat(struct2cell(snake_options_integer_valued));
  % snake_options_integer_valued = cellfun(@(x)snake_options_integer_valued.(x), snake_options_field_names);
  function [full_parameters] = usable_to_full_parameters(usable_parameters)
    full_parameters = reshape(cell2mat(struct2cell(initial_snake_options)), 1, []);
    % full_parameters(snake_options_to_use) = usable_parameters;
    % Disabled:
    full_parameters = usable_parameters;
  end
  function [usable_parameters] = full_to_usable_parameters(full_parameters)
    % usable_parameters = full_parameters(snake_options_to_use);
    % Disabled:
    usable_parameters = full_parameters;
  end
  function [usable_parameters] = snake_options_structure_to_vector(full_parameters)
    usable_parameters = reshape(cell2mat(struct2cell(full_parameters)), 1, []);
  end
  % final_snake_options = initial_snake_options;
    
  
  % Use evolutionary art approach:
  % Initialize:
  current_chromosome = struct2cell(initial_snake_options);
  snake_options_limits_cell = struct2cell(snake_options_limits);
  
  
  % This is originally from choose_bk_segmentation_histeq_usage_in_test_cells.m:
  function [evaluation_features] = compute_segmentation_evaluation_features(argument_structure)

    [preprocessing_structure] = compute_preprocessed_images(argument_structure);
    if length(fieldnames(preprocessing_structure)) == 0
      evaluation_features = [];
      return
    end
    synapse_locations_3d = preprocessing_structure.synapse_locations_3d;
    refined_synapse_locations_3d = preprocessing_structure.refined_synapse_locations_3d;
    dic_image = preprocessing_structure.dic_image;
    raw_image = preprocessing_structure.raw_image;
    x_range = preprocessing_structure.x_range;
    y_range = preprocessing_structure.y_range;
    z_range = preprocessing_structure.z_range;
    image_to_segment = preprocessing_structure.image_to_segment;
    coherent_image_to_segment2 = preprocessing_structure.coherent_image_to_segment2;

    % segmentation_structure = argument_structure.segmentation_structure;
    segmentation_structure = argument_structure.result;
    
    stretched_image = image_to_segment;
    
    segmentation_mesh = segmentation_structure.segmentation_mesh;
    % rotated_image = raw_image;
    % rotated_image = stretched_image;
    rotated_image = coherent_image_to_segment2;
    segmentation_image = voxelise(1:size(rotated_image, 2), 1:size(rotated_image, 1), 1:size(rotated_image, 3), segmentation_mesh);
    segmentation_image = permute(segmentation_image, [2, 1, 3]);
    rotated_segmentation = segmentation_image;
    
    % imshow([contrast_stretch(mean(rotated_image, 3)), contrast_stretch(mean(segmentation_image, 3))]), pause
    
    evaluation_features = [];
    % thresholded_rotated_image = im2bw(rotated_image);
    thresholded_rotated_image = rotated_image >= graythresh(rotated_image) * (max(rotated_image(:)) - min(rotated_image(:))) + min(rotated_image(:));
    rotated_segmentation_distance = bwdist(~rotated_segmentation);
    rotated_segmentation_distance_stretched = contrast_stretch(rotated_segmentation_distance);
    perimeter = rotated_segmentation_distance > 0 & rotated_segmentation_distance < 2;
    surface_area = sum(perimeter(:));
    volume = sum(rotated_segmentation(:));
    [gradient_x, gradient_y, gradient_z] = smooth_gradient(rotated_image);
    gradient_magnitude = sqrt(gradient_x.^2 + gradient_y.^2 + gradient_z.^2);
    [segmentation_gradient_x, segmentation_gradient_y, segmentation_gradient_z] = smooth_gradient(double(rotated_segmentation));
    segmentation_gradient_magnitude = sqrt(segmentation_gradient_x.^2 + segmentation_gradient_y.^2 + segmentation_gradient_z.^2);
    
    % Mean voxel intensity inside the cell divided by the mean outside:
    % proportion_nonfinite = mean(~isfinite(thresholded_rotated_image(:)))
    % mean(thresholded_rotated_image(rotated_segmentation))
    % mean(thresholded_rotated_image(~rotated_segmentation))
    evaluation_features = [evaluation_features, mean(thresholded_rotated_image(rotated_segmentation)) ./ mean(thresholded_rotated_image(~rotated_segmentation))];
    
    % Mean voxel intensity times (1 - distance from boundary / maximum distance) inside the cell:
    evaluation_features = [evaluation_features, mean(rotated_image(rotated_segmentation) .* (1 - rotated_segmentation_distance_stretched(rotated_segmentation)))];
    
    % Approximate surface area to volume ratio:
    evaluation_features = [evaluation_features, surface_area / volume];
    
    % Mean dot product of normalized gradients:
    % whos gradient_magnitude segmentation_gradient_magnitude
    usable_voxels = (gradient_magnitude > 0) & (segmentation_gradient_magnitude > 0);
    evaluation_features = [evaluation_features, mean((gradient_x(usable_voxels) .* segmentation_gradient_x(usable_voxels) + gradient_y(usable_voxels) .* segmentation_gradient_y(usable_voxels) + gradient_z(usable_voxels) .* segmentation_gradient_z(usable_voxels)) ./ (gradient_magnitude(usable_voxels) .* segmentation_gradient_magnitude(usable_voxels)))];
    % Mean powered/penalized inner product of gradients (feature 6):
    gradient_penalty = 1;
    % gradient_penalty = 2;
    evaluation_features = [evaluation_features, mean(((gradient_x(usable_voxels) .* segmentation_gradient_x(usable_voxels) + gradient_y(usable_voxels) .* segmentation_gradient_y(usable_voxels) + gradient_z(usable_voxels) .* segmentation_gradient_z(usable_voxels)) * -.5 + .5).^gradient_penalty)];
    
    % Mean powered/penalized inner product of gradients at vertices (feature 6):
    normals = PatchNormals3D(segmentation_mesh);
    vertices = segmentation_mesh.vertices;
    gradients = [...
      interp3(gradient_x, vertices(:, 1), vertices(:, 2), vertices(:, 3))...
      , interp3(gradient_y, vertices(:, 1), vertices(:, 2), vertices(:, 3))...
      , interp3(gradient_z, vertices(:, 1), vertices(:, 2), vertices(:, 3))...
      ];
    gradients(isnan(gradients)) = 0;
    evaluation_features = [evaluation_features, mean((sum(normals .* gradients, 2) * -.5 + .5).^gradient_penalty)];
    % evaluation_features
    
    % Mean intensity of interior divided by the same of the original ROI (for a better background sample):
    evaluation_features = [evaluation_features, mean(rotated_image(rotated_segmentation)) / mean(raw_image(:))];
    % Same, but using the median:
    % evaluation_features = [evaluation_features, median(rotated_image(rotated_segmentation)) / median(raw_image(:))];
    % Median is sometimes zero, so use normalized difference:
    evaluation_features = [evaluation_features, 1 - median(raw_image(:)) / median(rotated_image(rotated_segmentation))];
    
    evaluation_features(isnan(evaluation_features)) = 0;
    
    % keyboard
  end

  function [segmentation_quality] = minimizer_compute_segmentation_quality(current_parameters)
    current_parameters = usable_to_full_parameters(current_parameters);
    % Segment a random cell and automatically evaluate whether these parameters did better.
    number_parameter_sets = length(current_parameters) / length(snake_options_field_names);
    if number_parameter_sets == 1
      % current_parameters
      % snake_options = cell2struct(num2cell(current_parameters), snake_options_field_names);
      snake_options = cell2struct(num2cell(reshape(current_parameters, [], 1)), snake_options_field_names);
      snake_options.Verbose = false;
      % snake_options.Verbose = true;
      % gca
      snake_options.replicate_padding = true;
    else
      snake_options = {};
      for parameter_set_index = 1:number_parameter_sets
        current_parameters_subset = current_parameters((parameter_set_index - 1) * length(snake_options_field_names) + (1:length(snake_options_field_names)));
        % snake_options(parameter_set_index) = cell2struct(num2cell(reshape(current_parameters_subset, [], 1)), snake_options_field_names);
        snake_options{parameter_set_index} = cell2struct(num2cell(reshape(current_parameters_subset, [], 1)), snake_options_field_names);
        snake_options{parameter_set_index}.Verbose = false;
        % fprintf('>>>>>>>> HACK\n'), snake_options{parameter_set_index}.Verbose = true;
        % fprintf('>>>>>>>> HACK\n'), if parameter_set_index > 1, snake_options{parameter_set_index}.Verbose = true; end
        % gca
        snake_options{parameter_set_index}.replicate_padding = true;
      end
    end

    % Use the same random images:
    % Set random seed:
    % s = RandStream('mt19937ar', 'Seed', 672877);
    s = RandStream('mt19937ar', 'Seed', 70107);
    RandStream.setDefaultStream(s);

    options = base_options;
    options.dry_run = true;
    options.use_random_order = true;
    options.maximum_number_to_process = 10;
    % options.maximum_number_to_process = 20;
    % fprintf('>>>>>>>> HACK\n'), options.finished_compute_function = @segmentation_finished_compute_function;
    if fmin_use_single_image && optimization_running
      % fprintf('>>>>>>>> HACK was set to 1\n')
      options.maximum_number_to_process = 1;
    end
    % options.verbose = true;
    options.load_location = '/projects/cellorganizer/xruan/tcell_project/tcell/bhcho_bin_Tcell2/Results/Features_New_Initial_Set3_127/';
    options.save_location = '/projects/cellorganizer/xruan/tcell_project/tcell/bhcho_bin_Tcell2/Results/Features_New_Initial_Set3_127_refined_discard/';
    
    % segmentations = traverse_t_cell_original_images(...
      % @(v)segment_image(v, snake_options)...
      % , options);
    % Save memory:
    options.just_save_results = true;
    options.finished_compute_function = @(v)compute_segmentation_evaluation_features(v);
        
    if number_parameter_sets == 1
      [~, segmentation_qualities] = traverse_t_cell_original_images(...
        @(v)segment_image(v, snake_options)...
        , options);
    else
      % Use a two-stage segmentation method:
      % [~, segmentation_qualities] = traverse_t_cell_original_images(...
        % @(v)segment_image(setfield(v, 'previous_segmentation_mesh', getfield(segment_image(v, snake_options{1}), 'result')), snake_options{2})...
        % , options);
      [~, segmentation_qualities] = traverse_t_cell_original_images(...
        @(v)two_stage_segmentation(v)...
        , options);
    end
    
    % fprintf('segmentation_qualities:\n'), segmentation_qualities{:}
    
    % Evaluate segmentations using this objective function:
    segmentation_qualities = cell2mat(segmentation_qualities(~cellfun('isempty', segmentation_qualities)));
    % segmentation_quality = sum(segmentation_qualities(:, 7));
    % segmentation_quality = mean(segmentation_qualities(:, 7));
    segmentation_quality = mean(segmentation_qualities(:, 6));
    
    % Negate segmentation_quality for minimization:
    % current_parameters, segmentation_quality
    fprintf(['segmentation_quality = ', num2str(segmentation_quality), ': ', num2str(current_parameters, '%.4e, '), '\n']);
    segmentation_quality = -segmentation_quality;
  end
  
  function [stop] = diary_cycle(x, y, z)
    stop = false;
    diary off
    diary on
  end

  
  
  coarse_snake_options = struct(...
    'Gamma', 2e-1 ...
    , 'Sigma1', 1 ...
    , 'Wline', -5e0 ...
    , 'Wedge', 0 ...
    , 'Sigma2', 2 ...
    , 'Sigma3', 2 ...
    , 'Alpha', 3e-2 * 5 ...
    , 'Beta', 1e-1 * 1 ...
    , 'Delta', 1e-1 ...
    , 'Kappa', 1 * 1 ...
    , 'Iterations', 15 * 4 * 4 ...
    , 'GIterations', ceil(cell_radius * .5) * 0 ...
    , 'Lambda', .95 ...
    )
  fine_snake_options = struct(...
    'Gamma', 5e-1 * 10 ...
    , 'Sigma1', 1 ...
    , 'Wline', -5e0 ...
    , 'Wedge', 1e1 * 0 ...
    , 'Sigma2', 1 * sqrt(3)/2 * 0 + 1 ...
    , 'Sigma3', 1 * sqrt(3)/2 * 0 + 1 ...
    , 'Alpha', 1.5e-1 * 1e-0 ...
    , 'Beta', 1e-1 * 1e-0 ...
    , 'Delta', 1e-1 * 1e-0 ...
    , 'Kappa', 1 * 10 ...
    , 'Iterations', 15 * 4 * 4 * 4 / 4 ...
    , 'GIterations', ceil(cell_radius * .2 / .5) * 1 ...
    , 'Lambda', .75 ...
    )
  initial_fmincon_parameters = [snake_options_structure_to_vector(coarse_snake_options), snake_options_structure_to_vector(fine_snake_options)];
  
  % optimization_method = 'iga'
  % optimization_method = 'fmincon'
  % optimization_method = 'none', fprintf('>>>>>>>> HACK\n'), final_chromosome = cell2mat(struct2cell(initial_snake_options))';
  % optimization_method = 'none', fprintf('>>>>>>>> HACK\n'), final_chromosome = getfield(load([base_filename, 'r234pre_keep', filesep, 'iga_iter00005_indiv00001.mat']), 'individuals'); final_chromosome = final_chromosome(1, :);
  optimization_method = 'none', fprintf('>>>>>>>> HACK, using initial_fmincon_parameters without optimization\n'), final_chromosome = initial_fmincon_parameters;
  % optimization_method = 'none'
  

  
  switch optimization_method
    case 'fmincon'
      save_filename = [base_filename, optimization_method, '_full_output.mat'];
      if exist(save_filename, 'file')
        load(save_filename)
      else
        a = load([base_filename, 'r234pre_keep', filesep, 'iga_iter00004_indiv00016.mat']);
        bounds = cell2mat(snake_options_limits_cell);
        
        % Ignore parameters that will not change:
        % initial_parameters = initial_parameters(1, snake_options_to_use);
        % initial_parameters = full_to_usable_parameters(initial_parameters);
        % bounds
        % bounds = bounds(1, snake_options_to_use);
        % bounds = bounds(snake_options_to_use, 1);
        % Disabled:
        % bounds = bounds(snake_options_to_use, :);

        bounds = repmat(bounds, 2, 1);
        initial_parameters = initial_fmincon_parameters;
        
        % fprintf('>>>>>>>> HACK\n'), optimization_running = true, minimizer_compute_segmentation_quality(initial_parameters), optimization_running = false
        
        if matlabpool('size') > 0
          matlabpool('close');
        end
        matlabpool('open', 8 + feature('numCores') * 0);
        optimization_running = true
        try
          typical_x = (1e-3 / sqrt(eps)) ./ mean(bounds, 2) .* ~[snake_options_integer_valued; snake_options_integer_valued] + [snake_options_integer_valued; snake_options_integer_valued] ./ sqrt(eps);
          fmincon_options = optimset('Algorithm', 'interior-point', 'LargeScale', 'off', 'TypicalX', typical_x, 'FinDiffType', 'forward', 'DiffMinChange', 1e-3, 'DiffMaxChange', 2, 'TolFun', 1e-8, 'UseParallel', 'always', 'Display', 'iter', 'OutputFcn', @diary_cycle);
          [final_chromosome, final_value, exit_flag, output, lambda, grad, hessian] = fmincon(@minimizer_compute_segmentation_quality, initial_parameters, [], [], [], [], bounds(:, 1), bounds(:, 2), [], fmincon_options);
        catch err
          if matlabpool('size') > 0
            matlabpool('close');
          end
          rethrow(err);
        end
        if matlabpool('size') > 0
          matlabpool('close');
        end
        optimization_running = false
        final_chromosome = usable_to_full_parameters(final_chromosome);

        save(save_filename, 'initial_parameters', 'bounds', 'snake_options_field_names', 'final_chromosome', 'final_value', 'exit_flag', 'output', 'lambda', 'grad', 'hessian', 'fmincon_options')
      end
  end
  final_snake_options = arrayfun(@(x)cell2struct(num2cell(final_chromosome((x - 1) * length(snake_options_field_names) + (1:length(snake_options_field_names))))', snake_options_field_names), 1:round(length(final_chromosome) / length(snake_options_field_names)), 'UniformOutput', false);
  final_snake_options

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Run the segmentations or generate images for visual inspection (just_run_finished_compute_function = true):

  options = base_options;
  options.save_location = snake_segmentation_location;
  % Prevent OOM error:
  options.just_save_results = true;

  % Turn off for actual computation, turn on after everything's finished/running on head node:
  if just_run_finished_compute_function
    options.finished_compute_function = @segmentation_finished_compute_function;
    options.just_run_finished_compute_function = true;
    % options.use_random_order = true;
    % % Set random seed:
    % s = RandStream('mt19937ar', 'Seed', 7129);
  end
  
  % skip_segmentation = false;
  skip_segmentation = true;
  if skip_segmentation
    fprintf('>>>> HACK, segmentation turned off\n')
  else
  traverse_t_cell_original_images(...
    @(v)two_stage_segmentation(v, final_snake_options)...
    , options);
  end

  % error('Implementation yet unfinished below this line!')

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Save the region corresponding to each segmentation for later processing's convenience:

  function [regions_structure] = save_regions(argument_structure)
    % argument_structure
    if size(argument_structure.frame_synapse_locations, 1) == 0
      segmentation_structure = struct();
      return
    end

    [preprocessing_structure] = compute_preprocessed_images(argument_structure);
    if length(fieldnames(preprocessing_structure)) == 0
      segmentation_structure = struct();
      return
    end

    number_synapses = preprocessing_structure.number_synapses;
    synapse_locations_3d = preprocessing_structure.synapse_locations_3d;
    refined_synapse_locations_3d = preprocessing_structure.refined_synapse_locations_3d;
    dic_image = preprocessing_structure.dic_image;
    x_range = preprocessing_structure.x_range;
    y_range = preprocessing_structure.y_range;
    z_range = preprocessing_structure.z_range;
    image_to_segment = preprocessing_structure.image_to_segment;
    coherent_image_to_segment2 = preprocessing_structure.coherent_image_to_segment2;
    coarse_white_line_image_to_segment = preprocessing_structure.coarse_white_line_image_to_segment ;
    fine_white_line_image_to_segment = preprocessing_structure.fine_white_line_image_to_segment;
    
    regions_structure = struct();
    regions_structure.number_synapses = preprocessing_structure.number_synapses;
    regions_structure.synapse_locations_3d = preprocessing_structure.synapse_locations_3d;
    regions_structure.refined_synapse_locations_3d = preprocessing_structure.refined_synapse_locations_3d;
    regions_structure.dic_image = preprocessing_structure.dic_image;
    regions_structure.x_range = preprocessing_structure.x_range;
    regions_structure.y_range = preprocessing_structure.y_range;
    regions_structure.z_range = preprocessing_structure.z_range;
    regions_structure.image_to_segment = preprocessing_structure.image_to_segment;
    regions_structure.raw_image = preprocessing_structure.raw_image;
  end

  options = base_options;
  options.save_location = snake_region_location;
  % Prevent OOM error:
  options.just_save_results = true;
  % skip_region_saving = false;
  skip_region_saving = true;
  if skip_region_saving
    fprintf('>>>> HACK, region saving turned off\n')
  else
    traverse_t_cell_original_images(@(v)save_regions(v), options);
  end
  
  
  
  
  function [evaluation_image_filename] = get_evaluation_image_filename(argument_structure, filenames, synapse_index)
  % Get the filename for the PNG given in the evaluation CSV file:
    original_argument_structure = argument_structure;
    segmentation_structure = load(strrep(argument_structure.filename, snake_region_location, snake_segmentation_location));
    segmentation_structure = segmentation_structure.current_result;

    argument_structure = segmentation_structure.argument_structure;
    
    % Compute the image's filename and see if it is in evaluation_filenames:
    fixed_relative_path = regexprep(argument_structure.relative_path, '/+', '/');
    fixed_relative_path = fixed_relative_path(find(fixed_relative_path == '/', 1, 'first') + 1:find(fixed_relative_path == '/', 1, 'last') - 1);
    bk_filename = [strrep(fixed_relative_path, '/', '_'), num2str(argument_structure.frame_synapse_cell_ids(synapse_index), '_%d'), num2str(argument_structure.frame_synapse_locations(synapse_index, 3), '_%d'), num2str(argument_structure.frame_index, '_%d')];
    evaluation_image_filename = [bk_filename, '_original_refined', num2str(synapse_index, '_synapse%03d'), '.png'];
  end


  function [segmentation_features] = compute_mesh_segmentation_evaluation_features(argument_structure, filenames, labels)
  % Compute features to distinguish good and bad segmentations as labeled in the evaluation CSV file:
    
    % argument_structure, keyboard
    original_argument_structure = argument_structure;
    region_structure = argument_structure.input_structure.current_result;
    segmentation_structure = load(strrep(argument_structure.filename, snake_region_location, snake_segmentation_location));
    segmentation_structure = segmentation_structure.current_result;
    % region_structure, segmentation_structure, keyboard
    
    argument_structure = segmentation_structure.argument_structure;
    
    % Store everything as a cell array, then cell2mat at the end so we don't have to hard-code column indices:
    segmentation_features = {};
    
    for synapse_index = 1:segmentation_structure.number_synapses
      if nargin >= 2
        % Compute the image's filename and see if it is in evaluation_filenames:
        [segmentation_finished_compute_function_image_filename] = get_evaluation_image_filename(original_argument_structure, filenames, synapse_index);
        
        row_index = find(strcmp(segmentation_finished_compute_function_image_filename, filenames));
        
        if isempty(row_index)
          continue
        end
      end
        
      approximate_synapse_location_2d = argument_structure.frame_synapse_locations(synapse_index, 1:2);
      x_range = segmentation_structure.x_range{synapse_index};
      y_range = segmentation_structure.y_range{synapse_index};
      z_range = segmentation_structure.z_range{synapse_index};
      approximate_synapse_location_2d = approximate_synapse_location_2d - [x_range(1), y_range(1)];
      
      current_mesh = segmentation_structure.segmentation_mesh{synapse_index};
      current_region = region_structure.raw_image{synapse_index};
      
      current_x_range = 1:length(x_range);
      current_y_range = 1:length(y_range); 
      segmentation_image = voxelise(current_x_range, current_y_range, z_range, current_mesh); 
      segmentation_image = permute(segmentation_image, [2, 1, 3]); 
      
      show_debug_plot = false;
      % show_debug_plot = true;
      if show_debug_plot
        if nargin >= 3
          fprintf('True label is %d\n', labels(row_index))
        end
      
        current_region_contrast_stretched = contrast_stretch(current_region);
        imshow(reshape_2d([current_region_contrast_stretched; current_region_contrast_stretched .* ~segmentation_image], 4));
        pause
      end
      
      % Preparation for feature computation:
      
      vertex_intensities = interp3(current_region, current_mesh.vertices(:, 1), current_mesh.vertices(:, 2), current_mesh.vertices(:, 3));

      curvature_smoothing = 2;
      [minimum_curvature_direction, maximum_curvature_direction, minimum_curvature, maximum_curvature, mean_curvature, gaussian_curvature, normal] = compute_curvature(current_mesh.vertices, current_mesh.faces, struct('curvature_smoothing', curvature_smoothing, 'verb', false));
      minimum_curvature_direction = minimum_curvature_direction';
      maximum_curvature_direction = maximum_curvature_direction';
      normal = normal';
      
      neighbor_distances = build_euclidean_weight_matrix(triangulation2adjacency(current_mesh.faces), current_mesh.vertices, 0);
      neighbor_distances = neighbor_distances(find(neighbor_distances));

      [gradient_x, gradient_y, gradient_z] = smooth_gradient(current_region, 'scharr5', false, true);
      gradient_magnitude = sqrt(gradient_x.^2 + gradient_y.^2 + gradient_z.^2);
      vertex_gradient_magnitudes = interp3(gradient_magnitude, current_mesh.vertices(:, 1), current_mesh.vertices(:, 2), current_mesh.vertices(:, 3));
      vertex_gradients = [interp3(gradient_x, current_mesh.vertices(:, 1), current_mesh.vertices(:, 2), current_mesh.vertices(:, 3)), interp3(gradient_y, current_mesh.vertices(:, 1), current_mesh.vertices(:, 2), current_mesh.vertices(:, 3)), interp3(gradient_z, current_mesh.vertices(:, 1), current_mesh.vertices(:, 2), current_mesh.vertices(:, 3))];
      vertex_gradients(vertex_gradient_magnitudes > 0, :) = vertex_gradients(vertex_gradient_magnitudes > 0, :) ./ repmat(vertex_gradient_magnitudes(vertex_gradient_magnitudes > 0, :), 1, 3);
      
      gradient_scalar_products = sum(normal .* vertex_gradients, 2);
      
      [x, y, z] = meshgrid(current_x_range, current_y_range, z_range);
      % segmentation_convex_hull_image = ~isnan(tsearchn(current_mesh.vertices, delaunayn(current_mesh.vertices), [x(:), y(:), z(:)])); 
      % segmentation_convex_hull = convhulln(current_mesh.vertices);
      segmentation_convex_hull_image = reshape(inhull([x(:), y(:), z(:)], current_mesh.vertices), size(current_region)); 
      
      % interesting_percentiles = 0:25:100;
      interesting_percentiles = 0:10:100;
      
      % Feature computation:
      
      % Intensity distribution at vertices:
      segmentation_features{synapse_index, 1} = prctile(vertex_intensities, interesting_percentiles);
      % Curvature distribution at vertices:
      segmentation_features{synapse_index, 2} = [prctile(minimum_curvature, interesting_percentiles), prctile(maximum_curvature, interesting_percentiles)];
      % Length distribution of edges:
      segmentation_features{synapse_index, 3} = prctile(neighbor_distances, interesting_percentiles);
      % Gradient magnitude distribution at vertices:
      segmentation_features{synapse_index, 4} = prctile(vertex_gradient_magnitudes, interesting_percentiles);
      % Scalar product of image gradient direction and mesh normal distribution at vertices:
      segmentation_features{synapse_index, 5} = prctile(gradient_scalar_products, interesting_percentiles);
      % Segmentation image features:
      volume = sum(segmentation_image(:));
      solidity = volume / sum(segmentation_convex_hull_image(:));
      segmentation_features{synapse_index, 6} = [volume, solidity];
      % % Spherical harmonics of vertex components:
      % % :
      % % :
      % keyboard
    end
    
    segmentation_features = cell2mat(segmentation_features);
  end


  function [segmentation_labels] = compute_mesh_segmentation_evaluation_labels(argument_structure, filenames, labels)
  % Retrieve the labels indicating good and bad segmentations as in the evaluation CSV file:
    % argument_structure, keyboard
    original_argument_structure = argument_structure;
    region_structure = argument_structure.input_structure.current_result;
    segmentation_structure = load(strrep(argument_structure.filename, snake_region_location, snake_segmentation_location));
    segmentation_structure = segmentation_structure.current_result;
    % region_structure, segmentation_structure, keyboard
    
    argument_structure = segmentation_structure.argument_structure;
    
    % Store everything as a cell array, then cell2mat at the end so we don't have to hard-code column indices:
    segmentation_labels = {};
    
    for synapse_index = 1:segmentation_structure.number_synapses
      % Get this segmentation's image's filename and see if it is in evaluation_filenames:
      [segmentation_finished_compute_function_image_filename] = get_evaluation_image_filename(original_argument_structure, filenames, synapse_index);
      
      row_index = find(strcmp(segmentation_finished_compute_function_image_filename, filenames));
      
      if isempty(row_index)
        continue
      end
      
      if nargin >= 3
        segmentation_labels{synapse_index, 1} = labels(row_index);
      end
      
      use_debug = false;
      % use_debug = true;
      if use_debug
        fprintf('CSV file should say: %s, %d\n', segmentation_finished_compute_function_image_filename, labels(row_index));
      end
    end
    
    segmentation_labels = cell2mat(segmentation_labels);
  end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Filter out bad segmentations:
  
  
  evaluated_segmentation_features = [];
  evaluated_labels = [];
  evaluated_segmentation_features_filename = [base_filename, 'evaluated_segmentation_features.mat'];
  
  % This part has to run with use_test_cells = true:
  if use_test_cells
    % Training data for this SVM are labeled for test subset of cells:
    
    % Load manual labels:
    evaluation_csv_filename = 'segment_using_snakes_r266_evaluation.csv';
    a = importdata(evaluation_csv_filename, ',', 1);
    a.textdata = a.textdata(2:end, 1);
    a.textdata = strtrim(a.textdata);
    % a.textdata = strrep(a.textdata, '');
    evaluated_paths = a.textdata;
    evaluated_filenames = evaluated_paths;
    % evaluated_filenames
    evaluated_synapses = regexp(evaluated_filenames, '_original_refined_synapse([0-9]+)\.png', 'tokens');
    evaluated_synapses = cellfun(@(x)str2num(x{1}{1}), evaluated_synapses);
    evaluated_frames = regexp(evaluated_filenames, '_([0-9]+)_original_refined_synapse[0-9]+\.png', 'tokens');
    evaluated_frames = cellfun(@(x)str2num(x{1}{1}), evaluated_frames);
    evaluated_filenames = regexprep(evaluated_filenames, '(_([0-9]+)){4}_original_refined_synapse[0-9]+\.png', '');
    evaluated_filenames = strrep(evaluated_filenames, '_', '/');
    
    % evaluated_frames, evaluated_filenames, keyboard
    
    evaluated_file_ids = cellfun(@(x)t_cell_filename_to_file_id(x), a.textdata);
    number_evaluated_cells = length(evaluated_filenames)
    evaluated_cells_are_good = a.data;
    evaluated_labels = evaluated_cells_are_good + 1;
    
    % Get features just for labeled cells:
    
    % For traverse_directory_structure:
    options = struct();
    
    % % For traverse_t_cell_original_images:
    % options = base_options;
    
    options.dry_run = true;
    options.load_location = snake_region_location;
    options.save_location = snake_filtered_region_location;
    
    synapse_identifier_function = @(given_traverse_file_id)struct_ismember(rmfield(evaluated_file_ids, 'filename'), rmfield(t_cell_image_path_to_file_id(given_traverse_file_id, snake_region_location), 'filename'));
    
    options.compute_function = @(given_argument_structure)compute_mesh_segmentation_evaluation_features(given_argument_structure, evaluated_paths, evaluated_cells_are_good);
    
    options.finished_compute_function = @(given_argument_structure)compute_mesh_segmentation_evaluation_labels(given_argument_structure, evaluated_paths, evaluated_cells_are_good);
    
    options.should_compute_function = @(given_argument_structure)any(synapse_identifier_function(given_argument_structure.file_id));
    
    % % Prevent OOM error:
    % options.just_save_results = true;
    % Get labels in the same order as features instead of as they are in the CSV file:
    if ~exist(evaluated_segmentation_features_filename, 'file')
      [evaluated_segmentation_features, ~, evaluated_labels] = traverse_directory_structure(options);
      evaluated_segmentation_features = cell2mat(evaluated_segmentation_features);
      evaluated_labels = cell2mat(evaluated_labels);
      evaluated_labels = evaluated_labels + 1;
      save(evaluated_segmentation_features_filename, 'evaluated_segmentation_features', 'evaluated_labels')
    % else
      % load(evaluated_segmentation_features_filename)
    end
    % whos evaluated_segmentation_features evaluated_labels
    % keyboard
    
    % Get the original 2D manual synapse locations for these cells so we can identify them in the entire data set:
    warning('Need to get the original 2D manual synapse locations for these cells so we can identify them in the entire data set!!!!!!!!!!!!!!!!!!!!!')
  end
    

  if ~exist(evaluated_segmentation_features_filename, 'file')
    error('Run with use_test_cells = true first!')
  else
    load(evaluated_segmentation_features_filename)
  end
  
  % Train the SVM using cross-validated parameters:
  % fprintf('>>>> HACK, using noise features!\n'), evaluated_segmentation_features = rand(size(evaluated_segmentation_features));
  svm_info = train_svm_with_cross_validation(evaluated_segmentation_features, evaluated_labels, [], [base_filename, strrep(evaluation_csv_filename, '.csv', '_svm')], 'none', true);
  % keyboard
    


  function [cells_are_good] = evaluate_cells_using_svm(argument_structure)
  % Determine if each synapse is good:
    % number_synapses = argument_structure.input_structure.current_result.number_synapses;
    % for synapse_index = 1:number_synapses
    % corresponding_region = 
    original_features = compute_mesh_segmentation_evaluation_features(argument_structure);
    number_synapses = size(original_features, 1);
    standardized_features = (original_features - repmat(svm_info.classifier_info.training_mean, number_synapses, 1)) ./ repmat(svm_info.classifier_info.training_standard_deviation, number_synapses, 1);
    % predictions = svmpredict(zeros(number_synapses, 1), standardized_features, svm_info.classifier_info.final_model, '-b 1');
    predictions = svmpredict(zeros(number_synapses, 1), standardized_features, svm_info.classifier_info.final_model, '-b 1 -q');
    cells_are_good = predictions == 2;
  end
  

  
  function [segmentation_structure] = blank_synapses_using_svm(argument_structure)
    region_structure = load(strrep(argument_structure.filename, snake_segmentation_location, snake_region_location));
    region_structure = region_structure.current_result;
    segmentation_structure = argument_structure.input_structure.current_result;
    segmentation_structure.raw_image = region_structure.raw_image;
    
    feature_argument_structure = argument_structure;
    feature_argument_structure.input_structure.current_result = segmentation_structure;
    
    % keyboard
    % good_segmentations = evaluate_cells_using_svm(argument_structure);
    good_segmentations = evaluate_cells_using_svm(feature_argument_structure);

    % Blanking should be unnecessary due to saved segmentation_structure.good_segmentations, but this should cause errors down the line for anything that has not yet been fixed:
    segmentation_structure = argument_structure.input_structure.current_result;
    segmentation_structure.good_segmentations = good_segmentations;
    if ~all(good_segmentations)
      for bad_segmentation_index = reshape(find(~good_segmentations), 1, [])
        % keyboard
        segmentation_structure.segmentation_mesh{bad_segmentation_index} = [];
        segmentation_structure.synapse_locations_3d{bad_segmentation_index} = [];
        segmentation_structure.refined_synapse_locations_3d{bad_segmentation_index} = [];
        segmentation_structure.dic_image{bad_segmentation_index} = [];
        segmentation_structure.x_range{bad_segmentation_index} = [];
        segmentation_structure.y_range{bad_segmentation_index} = [];
        segmentation_structure.z_range{bad_segmentation_index} = [];
      end
    end
  end
  
    
  % Filter all cells based on SVM output:
  % options = base_options;
  options = struct();
  % options.dry_run = true;
  % options.use_random_order = true;
  options.load_location = snake_segmentation_location;
  % options.file_listing_function = @(given_load_location)dirr(given_load_location);
  % options.load_location = snake_region_location;
  options.save_location = snake_filtered_segmentation_location;
  % options.save_location = snake_filtered_region_location;
  % options.compute_function = @(given_compute_function_arguments)given_compute_function_arguments.input_structure;
  options.compute_function = @blank_synapses_using_svm;
  % error('This needs to use svmpredict and standardization to make a decision, put in a separate function!!!')
  % options.should_compute_function = @(given_compute_function_arguments)svmpredict(0, (compute_mesh_segmentation_evaluation_features(given_compute_function_arguments) - svm_info.classifier_info.training_mean) ./ svm_info.classifier_info.training_standard_deviation, svm_info.classifier_info.final_model, '-b 1') == 2;
  % options.should_compute_function = @evaluate_cells_using_svm;
  % Prevent OOM error:
  options.just_save_results = true;
  skip_segmentation_filtering = false;
  % skip_segmentation_filtering = true;
  if skip_segmentation_filtering
    fprintf('>>>> HACK, segmentation filtering turned off\n')
  else
    traverse_directory_structure(options);
  end
  % keyboard

  

  % % Haven't updated the below to account for segmentation filter yet (ideally just use snake_filtered_region_location instead of snake_region_location):
  % error('Implementation yet unfinished below this line!')

  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Align segmentations:
  
  
  function [aligned_segmentation_structure] = align_segmentation_structure(argument_structure)
    % Algorithmic alignment of individual frames.
    segmentation_structure = argument_structure.input_structure.current_result;
    argument_structure = segmentation_structure.argument_structure;

    aligned_segmentation_structure = segmentation_structure;
    aligned_segmentation_structure.raw_image_aligned = cell(segmentation_structure.number_synapses, 1);
    aligned_segmentation_structure.alignment_transforms = cell(segmentation_structure.number_synapses, 1);
    
    % save_debug_plots = false;
    save_debug_plots = true;
    
    align_options = struct();
    align_options.debug_plot = false;
    % align_options.debug_plot = true;
    show_image_debug_plot = false;
    % show_image_debug_plot = true;
    if save_debug_plots
      align_options.debug_plot = true;
    end

    % align_options.synapse_detection_method = 'curvature local minima';
    fprintf('>>>> NOTE, using align_options.synapse_detection_method = ''synapse direction''\n'), align_options.synapse_detection_method = 'synapse direction';

    
    x_range = segmentation_structure.x_range;
    y_range = segmentation_structure.y_range;
    z_range = segmentation_structure.z_range;
    
    align_options.images_to_transform = {};
    for cell_index = 1:segmentation_structure.number_synapses
      % Skip cells with bad segmentations according to the SVM:
      % x_range{cell_index}
      % x_range_isempty = isempty(x_range{cell_index})
      if isempty(x_range{cell_index})
        continue
      end
      
      fixed_relative_path = regexprep(argument_structure.relative_path, '/+', '/');
      fixed_relative_path = fixed_relative_path(find(fixed_relative_path == '/', 1, 'first') + 1:find(fixed_relative_path == '/', 1, 'last') - 1);
      bk_filename = [strrep(fixed_relative_path, '/', '_'), num2str(argument_structure.frame_synapse_cell_ids(cell_index), '_%d'), num2str(argument_structure.frame_synapse_locations(cell_index, 3), '_%d'), num2str(argument_structure.frame_index, '_%d')];

      % Get the corresponding image region:
      % error('Get the corresponding image region!')
      temp = load([snake_region_location, filesep, argument_structure.relative_path]);
      % keyboard
      cell_image_region = temp.current_result.raw_image{cell_index};

      segmentation_mesh = segmentation_structure.segmentation_mesh{cell_index};
      segmentation_result = segmentation_mesh;
      current_x_range = 1:length(x_range{cell_index});
      current_y_range = 1:length(y_range{cell_index});
      % Do not correct vertices to be in the space of the original image as we pass cell_image_region to align_synapse_in_segmentation_image:
      % segmentation_result.vertices(:, 1) = segmentation_result.vertices(:, 1) + x_range{cell_index}(1) - 1;
      % segmentation_result.vertices(:, 2) = segmentation_result.vertices(:, 2) + y_range{cell_index}(1) - 1;
      % current_x_range = current_x_range + x_range{cell_index}(1) - 1;
      % current_y_range = current_y_range + y_range{cell_index}(1) - 1;
      segmentation_image = voxelise(current_x_range, current_y_range, z_range{cell_index}, segmentation_result);
      segmentation_image = permute(segmentation_image, [2, 1, 3]);
      dilate_radius = 0;
      % dilate_radius = 2;
      intensity_mask = imdilate(segmentation_image, ones(ones(1, 3) * (1 + 2 * dilate_radius)));
      cell_image_region_masked = cell_image_region .* intensity_mask + min(cell_image_region(:)) .* ~intensity_mask;

      align_options.images_to_transform = {cell_image_region_masked};
      align_options.images_to_transform{end+1} = segmentation_image;
      
      % if align_options.debug_plot && ~save_debug_plots
      % if show_image_debug_plot && ~save_debug_plots
      if show_image_debug_plot
        % Add the manual synapse point to the segmentation and align that as well for display:
        segmentation_image_to_show = segmentation_image;
        
        % segmentation_structure_x_range = segmentation_structure.x_range
        % approximate_synapse_location_2d = argument_structure.frame_synapse_locations(1, 1:2);
        approximate_synapse_location_2d = argument_structure.frame_synapse_locations(cell_index, 1:2);
        % approximate_synapse_location_2d = approximate_synapse_location_2d - [segmentation_structure.x_range{1}(1), segmentation_structure.y_range{1}(1)];
        approximate_synapse_location_2d = approximate_synapse_location_2d - [segmentation_structure.x_range{cell_index}(1), segmentation_structure.y_range{cell_index}(1)];
    
        sx = round(approximate_synapse_location_2d(1));
        sy = round(approximate_synapse_location_2d(2));
        
        if save_debug_plots
          synapse_marker_radius = 1;
        else
          synapse_marker_radius = 2;
        end
        
        % % Skip if (sx, sy) is outside the cropped image:
        % Limit (sx, sy) to be inside the cropped image:
        sx = min(max(sx, synapse_marker_radius * 2 + 1), size(segmentation_image, 2) - synapse_marker_radius * 2);
        sy = min(max(sy, synapse_marker_radius * 2 + 1), size(segmentation_image, 1) - synapse_marker_radius * 2);
        
        segmentation_image_to_show(sy + (-synapse_marker_radius * 2:synapse_marker_radius * 2), sx + (-synapse_marker_radius * 2:synapse_marker_radius * 2), :) = false;
        segmentation_image_to_show(sy + (-synapse_marker_radius:synapse_marker_radius), sx + (-synapse_marker_radius:synapse_marker_radius), :) = true;
        align_options.images_to_transform{end+1} = segmentation_image_to_show;
      end

      % Find the synapse and align the segmentation mesh and several images along with it:
      
      % initial_mesh = segmentation_mesh;
      % Maybe this is why the image transformations aren't working:
      initial_mesh = segmentation_result;
      
      if align_options.debug_plot
        % clf
        close
      end
      
      % initial_mesh
      % [aligned_segmentation_structure.segmentation_mesh{cell_index}, aligned_segmentation_structure.raw_image_aligned{cell_index}, aligned_segmentation_structure.alignment_transforms{cell_index}] = align_synapse_in_segmentation_image(initial_mesh, approximate_synapse_location_2d, align_options);
      [aligned_segmentation_structure.segmentation_mesh{cell_index}, aligned_images, aligned_segmentation_structure.alignment_transforms{cell_index}] = align_synapse_in_segmentation_image(initial_mesh, approximate_synapse_location_2d([2, 1]), align_options);
      aligned_segmentation_structure.raw_image_aligned{cell_index} = aligned_images{1};
      aligned_segmentation_image = aligned_images{2};

      if align_options.debug_plot
        % colors = pmkmp(256, 'CubicL');
        colors = pmkmp(256, 'IsoL');
        colormap(colors);
      end
      
      if align_options.debug_plot && save_debug_plots
        set(gcf, 'Visible', 'off');
        image_filename = [base_filename, filesep, bk_filename, '_detected_synapse', num2str(cell_index, '%03d')];
        ppi = 90;
        % ppi = 150;
        % ppi = 300;
        export_fig(image_filename, '-opengl', '-png', '-a1', '-nocrop', ['-r', num2str(ppi)])
      end
      
      % When debugging, set the right colormap and pause:
      if align_options.debug_plot && ~save_debug_plots
      % if show_image_debug_plot && ~save_debug_plots
        pause
      end

      % if align_options.debug_plot && ~save_debug_plots
      if show_image_debug_plot
        colormap(gray);
        aligned_cell_image_region = aligned_images{1};
        aligned_segmentation_image_to_show = aligned_images{3};
        % whos cell_image_region segmentation_image_to_show aligned_cell_image_region aligned_segmentation_image_to_show
        % % Look for an error with dbstop:
        % discard = [contrast_stretch(cell_image_region); segmentation_image_to_show; contrast_stretch(aligned_cell_image_region); aligned_segmentation_image_to_show];

        % Plot the segmentation before and after standardization/alignment:
        % number_plot_rows = 5;
        % number_plot_rows = 3;
        number_plot_rows = 2;
        % clf, set(gcf, 'Visible', 'on'), imshow(reshape_contrast([segmentation_image_to_show; aligned_segmentation_image], 5)), pause
        % clf, set(gcf, 'Visible', 'on'), imshow(reshape_contrast([segmentation_image_to_show; aligned_images{3}], 5)), pause
        debug_image = reshape_contrast([contrast_stretch(cell_image_region); segmentation_image_to_show; contrast_stretch(aligned_cell_image_region); aligned_segmentation_image_to_show], 3);
        clf, imshow(debug_image)
        % colormap(pmkmp(256, 'CubicL'));
        colormap(gray);
        if ~save_debug_plots
          set(gcf, 'Visible', 'on'), 
          pause
        end
      end
      % if align_options.debug_plot
      if show_image_debug_plot && save_debug_plots
        set(gcf, 'Visible', 'off');
        image_filename = [base_filename, filesep, bk_filename, '_synapse_alignment', num2str(cell_index, '%03d')];
        % ppi = 90;
        % % ppi = 150;
        % % ppi = 300;
        % export_fig(image_filename, '-opengl', '-png', '-a1', '-nocrop', ['-r', num2str(ppi)])
        imwrite(debug_image, [image_filename, '.png'])
      end
    end

    % keyboard
    % pause
  end
  


  function [result_landmarks] = frame_landmark_function(given_segmentation_image, given_synapse_point, given_ranges)
    % Returns synapse and centroid:
    
    % Debugging to check if the synapse is reasonably placed:
    if given_synapse_point(1) > size(given_segmentation_image, 2) || given_synapse_point(2) > size(given_segmentation_image, 1) || any(given_synapse_point < 1)
      given_synapse_point = round([size(given_segmentation_image, 2), size(given_segmentation_image, 1)] ./ 2);
      warning('given_synapse_point is out of bounds of given_segmentation_image, setting to center of given_segmentation_image!')
    end
    % imshow(sum(given_segmentation_image, 3))
    % hold on
    % plot(given_synapse_point(1), given_synapse_point(2), 'rx')
    % % plot(given_synapse_point(1) - given_ranges(1, 1), given_synapse_point(2) - given_ranges(2, 1), 'b+')
    % hold off
    % pause

    
    % synapse_column = rasterize_mesh(given_mesh, struct('oversampling_scale', 4, 'cropping', {[given_synapse_point(1) + [-.5, .5]; given_synapse_point(2) + [-.5, .5]; given_ranges(3, :)]}));
    % whos given_segmentation_image, given_synapse_point
    % synapse_column = given_segmentation_image(given_synapse_point(1) + (-1:1), given_synapse_point(2) + (-1:1), given_ranges(3, :));
    % synapse_column = given_segmentation_image(given_synapse_point(1) + (-1:1), given_synapse_point(2) + (-1:1), given_ranges(3, 1):given_ranges(3, 2));
    % synapse_column_values = given_segmentation_image(given_synapse_point(1) + (-1:1), given_synapse_point(2) + (-1:1), given_ranges(3, 1):given_ranges(3, 2));
    synapse_column_values = given_segmentation_image(given_synapse_point(2) + (-1:1), given_synapse_point(1) + (-1:1), given_ranges(3, 1):given_ranges(3, 2));
    
    % synapse_column = synapse_column_values;
    % synapse_column = mean(synapse_column, 1);
    % synapse_column = mean(synapse_column, 2);
    % synapse_column = reshape(synapse_column, 1, []);
    % % synapse_z = mean(synapse_column .* (given_ranges(3, 1):given_ranges(3, 2)));
    % synapse_z = mean([given_ranges(3, 1), given_ranges(3, 2)]);
    % if sum(synapse_column) == 0
      % warning('Cell index current_cell_index has a nan Z coordinate, setting to half the image height!', current_cell_index)
    % else
      % synapse_column = synapse_column ./ sum(synapse_column);
      % synapse_z = synapse_column * (given_ranges(3, 1):given_ranges(3, 2))';
    % end

    % warning('frame_landmark_function changed for robustness, rerun alignments!', current_cell_index)
    
    [x, y, z] = meshgrid(1:size(given_segmentation_image, 2), 1:size(given_segmentation_image, 1), 1:size(given_segmentation_image, 3));
    % [x, y] = meshgrid(1:size(given_segmentation_image, 2), 1:size(given_segmentation_image, 1));
    horizontal_distances = sqrt((x - given_synapse_point(1)).^2 + (y - given_synapse_point(2)).^2);
    horizontal_weights = exp(-horizontal_distances);
    % horizontal_weights = horizontal_weights ./ sum(horizontal_weights(:));
    % horizontal_values = mean(given_segmentation_image, 3);
    segmentation_volume = sum(given_segmentation_image(:));
    % horizontal_values = horizontal_values ./ sum(horizontal_values(:));
    % horizontal_weights = horizontal_weights .* horizontal_values;
    % horizontal_weights = horizontal_weights ./ sum(horizontal_weights(:));
    % whos horizontal_weights given_segmentation_image
    % horizontal_weights = repmat(horizontal_weights, [1, 1, size(given_segmentation_image, 3)]) .* given_segmentation_image;
    horizontal_weights = horizontal_weights .* given_segmentation_image;
    horizontal_weights = horizontal_weights ./ sum(horizontal_weights(:));
    synapse_z = sum(z(:) .* horizontal_weights(:));
    
    
    centroid = sum([x(:) .* given_segmentation_image(:), y(:) .* given_segmentation_image(:), z(:) .* given_segmentation_image(:)], 1) ./ segmentation_volume;
    
    
    
    
    % imshow(reshape_contrast(rasterize_mesh(given_mesh, struct('oversampling_scale', 4)), 4))
    % keyboard
    % size(synapse_column), size(given_ranges(3, 1):given_ranges(3, 2))
    % result_landmarks = [given_synapse_point, mean(synapse_column .* (given_ranges(3, 1):given_ranges(3, 2)))];
    % whos given_synapse_point synapse_z
    result_landmarks = [given_synapse_point, synapse_z];
    result_landmarks = [result_landmarks; centroid];
    
    
    % % Debugging to see if this looks like the centroid:
    % imshow(reshape_contrast(given_segmentation_image, floor(sqrt(size(given_segmentation_image, 3))))), pause
    % plot(synapse_column', 'b-'), hold on, plot([1, 1] * result_landmarks(3), [0, 1], 'r:'), hold off, pause
    % plot(squeeze(sum(sum(horizontal_weights, 1), 2)), 'b-'), hold on, plot([1, 1] * result_landmarks(3), [0, 1], 'r:'), hold off, pause
    
  end



  
  use_individual_synapse_alignment = false
  % use_individual_synapse_alignment = true
  
  
  if use_individual_synapse_alignment
  
    % Find individual synapses and reorient the segmentation:
    
    % options = base_options;
    options = struct();
    options.dry_run = true;
    options.use_random_order = true;
    % options.load_location = snake_segmentation_location;
    options.load_location = snake_filtered_segmentation_location;
    options.save_location = snake_aligned_segmentation_location;
    options.compute_function = @align_segmentation_structure;
    % Prevent OOM error:
    options.just_save_results = true;
    skip_synapse_alignment = false;
    % skip_synapse_alignment = true;
    if skip_synapse_alignment
      fprintf('>>>> HACK, synapse alignment turned off\n')
    else
      traverse_directory_structure(options);
    end
  
  else
  
    % Align all segmentations of a time series:
    % error('Implementation yet unfinished below this line!')
    
    % Get a list of all synapse files:
    % synapse_file_list = dirr(base_options.synapse_location, 'name', '.+\.[^.]+');
    % synapse_file_list = dirr([base_options.synapse_location, filesep, '*.txt']);
    % synapse_file_list = dirr(base_options.synapse_location, filesep, '*.txt']);
    % {synapse_file_list.name}.'
    % [~, ~, synapse_file_list] = dirr([base_options.synapse_location, filesep, '*.txt'], 'name');
    [~, ~, synapse_file_list] = dirr(base_options.synapse_location, 'name', 'synapse.txt');
    % synapse_file_list'
    % keyboard
    
    % % Collect synapse info:
    number_synapse_files = length(synapse_file_list);
    
    show_debug_plot = false;
    % show_debug_plot = true;
    
    % Loop through all time series:
    for synapse_file_index = 1:number_synapse_files
      synapse_file = synapse_file_list{synapse_file_index};
      synapse_image_path = strrep(fileparts(synapse_file), base_options.synapse_location, base_options.image_location);
      synapse_image_path = regexprep(synapse_image_path, '/+', '/');
      synapse_image_path
      synapse_image_list = dir([synapse_image_path, filesep, 'T*C*Z.stk']);
      synapse_image_list = {synapse_image_list(:).name}.';
      % synapse_image_list
      synapse_locations = get_synapse_file_data(synapse_file_list{synapse_file_index}, true);
      % keyboard
      number_cells = max(cell2mat(reshape(synapse_locations, [], 1)), [], 1);
      number_cells = number_cells(4);
      cell_ids = cell2mat(reshape(synapse_locations, [], 1));
      cell_ids = cell_ids(:, 4);
      % number_cells, cell_ids, pause
      number_frames = length(synapse_locations);
      % We can store all the structures in memory (one cell is ~100 KB):
      cell_segmentation_tracks = cell(number_cells, number_frames);
      cell_synapse_point_tracks = cell(number_cells, number_frames);
      cell_range_tracks = cell(number_cells, number_frames);
      cell_transformation_tracks = repmat({eye(4)}, [number_cells, number_frames]);
      cell_transformation_affine_tracks = cell_transformation_tracks;
      
      if show_debug_plot
        % close
        clf
        % camlight
        % axis equal
        hold on
        % patch_colors = squeeze(hsv2rgb((0:number_cells - 1)' ./ number_cells, .9 * ones(number_cells, 1), .9 * ones(number_cells, 1)));
        % patch_colors = squeeze(hsv2rgb((0:number_frames - 1)' ./ number_frames, .9 * ones(number_frames, 1), .9 * ones(number_frames, 1)));
        % patch_colors = pmkmp(number_frames, 'IsoL');
        % colormap(pmkmp(256, 'IsoL'));
      end
      
      % Loop through all frames and assign them to one time series per cell:
      maximum_region_size = [0, 0, 0];
      for frame_index = 1:number_frames
        % Add segmentation for this frame to the appropriate cell:
        frame_synapse_locations = synapse_locations{frame_index};
        frame_image_filename = [strrep(fileparts(synapse_file), base_options.synapse_location, snake_filtered_segmentation_location), filesep, strrep(synapse_image_list{frame_index}, '.stk', '.mat')];
        if ~exist(frame_image_filename, 'file')
          continue
        end
        frame_segmentation_structure = load(frame_image_filename);
        frame_segmentation_structure = frame_segmentation_structure.current_result;
        % frame_segmentation_structure
        % keyboard
        % for frame_cell_index = 1:size(frame_synapse_locations, 1)
        for frame_cell_index = 1:min(size(frame_synapse_locations, 1), length(frame_segmentation_structure.x_range))
          current_cell_index = frame_synapse_locations(frame_cell_index, 4);
          % current_cell_index
          if isempty(frame_segmentation_structure.x_range{frame_cell_index})
            continue
          end
          maximum_region_size = max(maximum_region_size, [range(frame_segmentation_structure.x_range{frame_cell_index}), range(frame_segmentation_structure.y_range{frame_cell_index}), range(frame_segmentation_structure.z_range{frame_cell_index})]);
          cell_segmentation_tracks{current_cell_index, frame_index} = frame_segmentation_structure.segmentation_mesh{frame_cell_index};
          cell_synapse_point_tracks{current_cell_index, frame_index} = frame_synapse_locations(frame_cell_index, 1:2);
          cell_range_limits_tracks{current_cell_index, frame_index} = [frame_segmentation_structure.x_range{frame_cell_index}([1, end]); frame_segmentation_structure.y_range{frame_cell_index}([1, end]); frame_segmentation_structure.z_range{frame_cell_index}([1, end])];
          
          if show_debug_plot
            current_mesh_to_display = cell_segmentation_tracks{current_cell_index, frame_index};
            current_mesh_to_display_offset = [size(frame_segmentation_structure.dic_image{frame_cell_index}), 0] .* .5 .* [current_cell_index, 0, 0];
            current_mesh_to_display.vertices = current_mesh_to_display.vertices + repmat(current_mesh_to_display_offset, size(current_mesh_to_display.vertices, 1), 1);
            patch(current_mesh_to_display, 'CData', frame_index, 'CDataMapping', 'scaled', 'FaceColor', 'flat', 'EdgeColor', 'none');
          end

        end
        % pause
      end

      % rasterization_oversampling_scale = 1;
      rasterization_oversampling_scale = 2;
      % rasterization_oversampling_scale = 4;
      segmentation_rasterization_function = @(given_mesh)rasterize_mesh(given_mesh, struct('oversampling_scale', rasterization_oversampling_scale, 'cropping', {[ones(3, 1), maximum_region_size']}));
      % segmentation_landmark_function = @(given_mesh)rasterize_mesh(given_mesh, struct('oversampling_scale', rasterization_oversampling_scale, 'cropping', {[ones(3, 1), maximum_region_size']}));
      
      % Align each cell's time series separately (first run for multiple values of landmark_relative_cost, then decide on value that best respects both types of error):
      % TODO: convert this to an automatic random sampling of cells from different cell lines:
      debug_determine_best_landmark_relative_cost = false;
      % debug_determine_best_landmark_relative_cost = true;
      if debug_determine_best_landmark_relative_cost
        % landmark_relative_costs = 2.^(-12:1:8);
        % landmark_relative_costs = 2.^(-12:1:12);
        landmark_relative_costs = 2.^(-12:1:20);
      else
        % error('Must determine best cost before running with debug_determine_best_landmark_relative_cost false!')
        landmark_relative_costs = 2^-6;
      end
      final_transform_translation_norms = nan(1, length(landmark_relative_costs));
      final_transform_rotation_norms = nan(1, length(landmark_relative_costs));
      for cost_index = 1:length(landmark_relative_costs)
        landmark_relative_cost = landmark_relative_costs(cost_index);
        % cost_index, landmark_relative_cost
        fprintf('landmark_relative_cost = %e\n', landmark_relative_cost);
      
        for current_cell_index = 1:number_cells
          cell_frames = ~cellfun(@isempty, cell_segmentation_tracks(current_cell_index, :));
          number_segmented_frames = sum(cell_frames);
          
          cell_segmentations = cell_segmentation_tracks(current_cell_index, cell_frames);
          cell_segmentation_images = cell(1, number_segmented_frames);
          cell_synapse_points = cell_synapse_point_tracks(current_cell_index, cell_frames);
          cell_cropped_synapse_points_3d = cell(1, number_segmented_frames);
          cell_range_limits = cell_range_limits_tracks(current_cell_index, cell_frames);
          
          cell_aligned_segmentations = cell_segmentations;
          cell_aligned_segmentation_images = cell_segmentation_images;
          % cell_aligned_synapse_points = cell_synapse_points;
          cell_aligned_cropped_synapse_points_3d = cell_cropped_synapse_points_3d;

          % cell_filename_prefix = [base_filename, strrep(strrep(fileparts(synapse_file), [base_options.synapse_location, '/'], ''), '/', '_'), num2str(current_cell_index, '_cell%05d_')];
          cell_filename_prefix = [base_filename, strrep(strrep(fileparts(synapse_file), [base_options.synapse_location, '/'], ''), '/', '_'), num2str(landmark_relative_cost, '_landmark-weight%1.1e'), num2str(current_cell_index, '_cell%05d_')];
          cell_transforms_filename = [cell_filename_prefix, 'transforms.mat'];
          [~, cell_transforms_filename_without_directory] = fileparts(cell_transforms_filename);
          

          [can_start, final_name, final_exists] = chunk_start(base_filename, cell_transforms_filename_without_directory);
          
          if can_start
            % Run the alignment and save the returned parameters:
            
            % Rasterize each frame's mesh:
            for frame_index = 1:number_segmented_frames
              % cell_segmentation_images{frame_index} = segmentation_rasterization_function(frame_index);
              cell_segmentation_images{frame_index} = segmentation_rasterization_function(cell_segmentations{frame_index});
            end
            
            % Function to get 3D synapse and centroid points for each frame:
            cell_frame_landmark_function = @(image_index)frame_landmark_function(cell_segmentation_images{image_index}, cell_synapse_points{image_index} - [cell_range_limits{image_index}(1, 1), cell_range_limits{image_index}(2, 1)] + 1, cell_range_limits{image_index});

            % Produce each frame's 3D synapse in the cropped image:
            for frame_index = 1:number_segmented_frames
              cell_cropped_synapse_points_3d{frame_index} = cell_frame_landmark_function(frame_index);
              cell_cropped_synapse_points_3d{frame_index} = cell_cropped_synapse_points_3d{frame_index}(1, :);
            end

            align_options = struct;
            align_options.number_images = number_segmented_frames;
            align_options.image_function = @(image_index)cell_segmentation_images{image_index};
            align_options.landmark_function = cell_frame_landmark_function;
            align_options.landmark_relative_cost = landmark_relative_cost;
            align_options.maximum_function_evaluations = number_segmented_frames * 5 * maximum_region_size(1) * 5;
            align_options.return_final_errors = true;
            if debug_determine_best_landmark_relative_cost
              align_options.absolute_error_image_filename_prefix = cell_filename_prefix;
              align_options.debug = 1;
            end
            
            % Could write these by resaving segmentations after transformation at some point:
            cell_alignment_transforms = align_time_series(align_options);
            % save(cell_transforms_filename, 'cell_alignment_transforms')
            
            mean_image_error = cell_alignment_transforms.final_image_error / number_segmented_frames;
            mean_landmark_error = cell_alignment_transforms.final_landmark_error / number_segmented_frames;
            
            save(final_name, 'cell_alignment_transforms', 'number_segmented_frames', 'landmark_relative_cost', 'current_cell_index', 'synapse_file', 'mean_image_error', 'mean_landmark_error')
            
            chunk_finish(base_filename, cell_transforms_filename_without_directory);
            
          % elseif final_exists
            % % Write image regions and segmentations to a new directory after transformation:
          
            % a = load(final_name);
            % cell_alignment_transforms = a.cell_alignment_transforms
            
            % for frame_index = 1:number_segmented_frames
              % % Align segmentation:
              % cell_aligned_segmentations{frame_index} = cell_segmentations{frame_index};
              % frame_rotations = sum(cell_alignment_transforms.rotations(1:frame_index - 1, :), 1);
              % frame_translations = sum(cell_alignment_transforms.translations(1:frame_index - 1, :), 1);
              % cell_aligned_segmentations{frame_index}.vertices = rigidly_transform_points(cell_aligned_segmentations{frame_index}.vertices, size(cell_segmentation_images{frame_index}), frame_rotations, frame_translations);
              % cell_aligned_segmentation_images{frame_index} = segmentation_rasterization_function(cell_aligned_segmentations{frame_index});
              
              % % Align synapse:
              % cell_aligned_cropped_synapse_points_3d{frame_index} = rigidly_transform_points(cell_cropped_synapse_points_3d{frame_index}, size(cell_segmentation_images{frame_index}), frame_rotations, frame_translations);
            % end
            
            % error('Write image regions and segmentations to a new directory after transformation')
            
          end
            
          % keyboard
        end
      end
      
      
      if debug_determine_best_landmark_relative_cost
        % plot landmark and image errors:
        addpath(genpath('boundedline'))
        set(0, 'DefaultAxesFontName', 'Helvetica')
        set(0, 'DefaultAxesFontSize', 16)
        set(0, 'DefaultLineLineWidth', 3)
        rd = 'segment_using_snakes/';
        %landmark_relative_costs = 2.^(-12:1:8);
        %landmark_relative_costs = 2.^(-12:1:9);
        landmark_relative_costs = 2.^(-12:1:20);
        image_errors = cell(length(landmark_relative_costs), 1);
        landmark_errors = cell(length(landmark_relative_costs), 1);
        maximum_index = length(landmark_relative_costs);
        %maximum_index = length(landmark_relative_costs) - 10;
        for cost_index = 1:maximum_index
            landmark_relative_cost = landmark_relative_costs(cost_index);
            fl = dir([rd, 'ARP3_050409 ARP3_Run 1', num2str(landmark_relative_cost, '_landmark-weight%1.1e'), '*_transforms.mat']);
            fprintf('cost %e: %d files\n', landmark_relative_cost, length(fl))
            for file_index = 1:length(fl)
                a = load([rd, fl(file_index).name]);
                image_errors{cost_index} = [image_errors{cost_index}, a.cell_alignment_transforms.final_image_error];
                landmark_errors{cost_index} = [landmark_errors{cost_index}, a.cell_alignment_transforms.final_landmark_error];
            end
        end
        mean_image_errors = cellfun(@(x)mean(x), image_errors)';
        mean_landmark_errors = cellfun(@(x)mean(x), landmark_errors)';
        interval_function = @(x)mean(x) + [-1, 1] .* norminv(1 - .05/2) .* std(x) ./ sqrt(length(x));
        image_error_intervals = cell2mat(cellfun(interval_function, image_errors, 'UniformOutput', false))';
        landmark_error_intervals = cell2mat(cellfun(interval_function, landmark_errors, 'UniformOutput', false))';
        good_indices = ~isnan(mean_image_errors);

        close, figure('Position', [1, 1, 8, 6] .* get(0, 'ScreenPixelsPerInch'))
        % subplot(1, 2, 1)
        hi1 = loglog(landmark_relative_costs(1:maximum_index), mean_image_errors(1:maximum_index), 'r:');
        hold on
        hl1 = loglog(landmark_relative_costs(1:maximum_index), mean_landmark_errors(1:maximum_index), 'b:');
        hi2 = loglog(landmark_relative_costs(1:maximum_index), image_error_intervals(:, 1:maximum_index)', 'r:')
        hl2 = loglog(landmark_relative_costs(1:maximum_index), landmark_error_intervals(:, 1:maximum_index)', 'b:')
        hold off

        axis tight
        % legend('Mean image error (w/95% CI)', 'Mean landmark error (w/95% CI)', 'Location', 'best')
        xlabel('Landmark weight')
        ylabel('Error')
        
        smoothing_kernel = sum(sum(gaussian_kernel(2), 3), 1);
        smoothing_kernel = smoothing_kernel ./ sum(smoothing_kernel(:));
        smoothed_mean_image_errors = convnfft_fast(mean_image_errors, smoothing_kernel, 'replicate')
        smoothed_mean_landmark_errors = convnfft_fast(mean_landmark_errors, smoothing_kernel, 'replicate')
        % smoothed_mean_image_error_slopes = diff(smoothed_mean_image_errors) ./ smoothed_mean_image_errors
        % smoothed_landmark_image_error_slopes = diff(smoothed_mean_landmark_errors) ./ smoothed_mean_landmark_errors
        smoothed_mean_image_error_slopes = convnfft_fast(smoothed_mean_image_errors, [1, 0, -1] / 2, 'replicate') ./ smoothed_mean_image_errors
        smoothed_landmark_image_error_slopes = convnfft_fast(smoothed_mean_landmark_errors, [1, 0, -1] / 2, 'replicate') ./ smoothed_mean_landmark_errors
        smoothed_relative_slope_sums = smoothed_mean_image_error_slopes + smoothed_landmark_image_error_slopes;
        [~, maximum_sum_index] = max(smoothed_relative_slope_sums)
        optimal_landmark_relative_cost = landmark_relative_costs(maximum_sum_index)
        
        % Plot optimal value:
        % subplot(1, 2, 2)
        hold on
        hi3 = loglog(landmark_relative_costs(1:maximum_index), smoothed_mean_image_errors(1:maximum_index), 'r-');
        hl3 = loglog(landmark_relative_costs(1:maximum_index), smoothed_mean_landmark_errors(1:maximum_index), 'b-');
        ho = loglog(optimal_landmark_relative_cost * [1, 1], ylim, 'g--');
        legend([hi1, hl1, hi3, hl3, ho], {'Mean image error (w/95% CI)', 'Mean landmark error (w/95% CI)', 'Smoothed mean image error', 'Smoothed mean landmark error', 'Optimal landmark relative cost'}, 'Location', 'best', 'Interpreter', 'none')
        hold off

        bk_filename = strrep(strrep(fileparts(synapse_file), [base_options.synapse_location, '/'], ''), '/', '_');
        image_filename = [base_filename, bk_filename, '_landmark_relative_cost_effects.png'];
        ppi = 90;
        % ppi = 150;
        % ppi = 300;
        export_fig(image_filename, '-opengl', '-png', '-a1', '-nocrop', ['-r', num2str(ppi)])
        
        return
      end
      

      % % Create spherical coordinate systems based on segmentation geometry and alignment to standardize cells:
      
      
      % Smooth the orientations of frames as found by fit_plane_to_synapse_in_mesh and store the final transformations and segmentations for morphing:
      for current_cell_index = 1:number_cells
        % current_title = sprintf('Time series ''%s'', cell %d', synapse_image_path, current_cell_index);
        current_title = sprintf('Time series ''%s'', cell %d', strrep(synapse_image_path, base_options.image_location, ''), current_cell_index);
        fprintf([current_title, '\n'])
        cell_frames_logical = ~cellfun(@isempty, cell_segmentation_tracks(current_cell_index, :));
        cell_frames = find(cell_frames_logical);
        cell_segmentations = cell_segmentation_tracks(current_cell_index, cell_frames_logical);
        cell_raw_images = cell(size(cell_segmentations));
        cell_synapse_points = cell_synapse_point_tracks(current_cell_index, cell_frames_logical);
        cell_range_limits = cell_range_limits_tracks(current_cell_index, cell_frames_logical);
        
        % % cell_filename_prefix = [base_filename, strrep(strrep(fileparts(synapse_file), base_options.synapse_location, ''), '/', '_'), num2str(current_cell_index, '_cell%05d_')];
        % cell_filename_prefix = [base_filename, strrep(strrep(fileparts(synapse_file), [base_options.synapse_location, '/'], ''), '/', '_'), num2str(current_cell_index, '_cell%05d_')];
        % cell_transforms_filename = [cell_filename_prefix, 'transforms.mat'];
        % [~, cell_transforms_filename_without_directory] = fileparts(cell_transforms_filename);
        
        [can_start, final_name, final_exists] = chunk_start(base_filename, cell_transforms_filename_without_directory);
        % if ~final_exists
          % chunk_finish(base_filename, cell_transforms_filename_without_directory);
        % else
          number_segmented_frames = length(cell_segmentations);
          cell_alignment_transforms = [];
          a = load(final_name)
          load(final_name, 'cell_alignment_transforms')
          if ~isfield(a, 'cell_alignment_transforms')
            % keyboard
            % pause
            continue
          end
          
          segmentation_images = cell(1, number_segmented_frames);
          cell_segmentations_aligned = cell(size(cell_segmentations));
          cell_synapse_points_aligned = zeros(size(cell_segmentations, 1), 3);
          for frame_index = 1:number_segmented_frames
            current_segmentation = cell_segmentations{frame_index};
            % segmentation_images{frame_index} = segmentation_rasterization_function(current_segmentation);
            
            segmentation_image_without_alignment = segmentation_rasterization_function(current_segmentation);
            current_synapse_point = frame_landmark_function(segmentation_image_without_alignment, cell_synapse_points{frame_index} - [cell_range_limits{frame_index}(1, 1), cell_range_limits{frame_index}(2, 1)] + 1, cell_range_limits{frame_index});
            % Centroid may be second landmark:
            % current_synapse_point = current_synapse_point(:, 1);
            current_synapse_point = current_synapse_point(1, :);
            
            % Transform according to alignment:
            total_rotation = sum(cell_alignment_transforms.rotations, 1);
            total_translation = sum(cell_alignment_transforms.translations, 1);
            % segmentation_images{frame_index} = rigidly_transform_image(segmentation_images{frame_index}, total_rotation, total_translation);
            point_transform_function = @(given_points)rigidly_transform_points(given_points, maximum_region_size, total_rotation, total_translation);
            
            % Accumulate transformations with respect to cropped original image:
            [~, transform, transform_affine] = point_transform_function(zeros(1, 3));
            cell_transformation_tracks{current_cell_index, frame_index} = transform * cell_transformation_tracks{current_cell_index, frame_index};
            cell_transformation_affine_tracks{current_cell_index, frame_index} = transform_affine * cell_transformation_affine_tracks{current_cell_index, frame_index};
            
            current_segmentation.vertices = point_transform_function(current_segmentation.vertices);
            current_synapse_point = point_transform_function(current_synapse_point);
            % cell_synapse_points_aligned{frame_index} = current_synapse_point;
            cell_synapse_points_aligned(frame_index, :) = current_synapse_point;
            segmentation_images{frame_index} = segmentation_rasterization_function(current_segmentation);
            cell_segmentations_aligned{frame_index} = current_segmentation;
          end
          
          % Smooth synapse position across time (note that frames are not evenly spaced as yet!):
          % cell_synapse_points_aligned_smoothed = convn(cell_synapse_points_aligned_smoothed, f', 'same');
          % temporal_smoothing_sigma = 1;
          % temporal_smoothing_sigma = 2;
          % Effectively turn off the per-frame effects of synapse orientation and position detection to see how the rigid alignment performed:
          temporal_smoothing_sigma = number_segmented_frames * 2;
          temporal_smoothing_filter = normpdf(-ceil(temporal_smoothing_sigma * 3):ceil(temporal_smoothing_sigma * 3), 0, temporal_smoothing_sigma);
          % This function smooths along columns:
          temporal_smoothing_function = @(given_sequence)convnfft_fast(given_sequence, temporal_smoothing_filter', 'replicate');
          
          cell_synapse_points_aligned_smoothed = temporal_smoothing_function(cell_synapse_points_aligned);
          centroids = zeros(number_segmented_frames, 3);
          % mean_centroid = [0, 0, 0];
          % maximum_region_size
          % [x, y, z] = meshgrid(1:maximum_region_size(1), 1:maximum_region_size(2), 1:maximum_region_size(3));
          % This is one larger in all dimensions for some reason... need to fix rasterize_mesh at some point:
          [x, y, z] = meshgrid(1:size(segmentation_images{1}, 2), 1:size(segmentation_images{1}, 1), 1:size(segmentation_images{1}, 3));
          for frame_index = 1:number_segmented_frames
            current_segmentation_image = segmentation_images{frame_index};
            % whos current_segmentation_image
            current_segmentation_volume = sum(current_segmentation_image(:));
            centroids(frame_index, :) = sum([x(:) .* current_segmentation_image(:), y(:) .* current_segmentation_image(:), z(:) .* current_segmentation_image(:)], 1) ./ current_segmentation_volume;
          end
          mean_centroid = mean(centroids, 1);
          centroids_smoothed = temporal_smoothing_function(centroids);
          
          % Compute synapse plane parameters and smooth them temporally:
          synapse_detection_results = cell(number_segmented_frames, 1);
          for frame_index = 1:number_segmented_frames
            current_segmentation = cell_segmentations_aligned{frame_index};
            current_segmentation_image = segmentation_images{frame_index};
            current_synapse_point = cell_synapse_points_aligned(frame_index, :);
            % align_synapse_in_segmentation_image(initial_mesh, approximate_synapse_location_2d([2, 1]), align_options);
            synapse_detection_results{frame_index} = fit_plane_to_synapse_in_mesh(current_segmentation, current_synapse_point, struct('synapse_detection_method', 'synapse direction'));
          end
          synapse_detection_results = cell2mat(synapse_detection_results);
          synapse_centers = cat(1, synapse_detection_results(:).synapse_center);
          synapse_normals = cat(1, synapse_detection_results(:).synapse_normal);
          synapse_radii = cat(1, synapse_detection_results(:).synapse_radius);
          synapse_centers_smoothed = temporal_smoothing_function(synapse_centers);
          synapse_normals_smoothed = temporal_smoothing_function(synapse_normals);
          synapse_normals_smoothed = synapse_normals_smoothed ./ repmat(sqrt(sum(synapse_normals_smoothed.^2, 2)), 1, 3);
          synapse_radii_smoothed = temporal_smoothing_function(synapse_radii);
          
          % synapse_centers, synapse_normals, synapse_radii
          % synapse_centers_smoothed, synapse_normals_smoothed, synapse_radii_smoothed
          
          template_parameters = get_default_template_options();
          template_image = get_template(template_parameters);
          [x, y, z] = meshgrid(1:size(template_image, 2), 1:size(template_image, 1), 1:size(template_image, 3));
          template_image_centroid = sum([x(:) .* template_image(:), y(:) .* template_image(:), z(:) .* template_image(:)], 1) ./ sum(template_image(:));

          % Replace cell_segmentations_aligned's entries with the segmentations transformed to align with the template using the smoothed parameter results from immediately above:
          for segmented_frame_index = 1:number_segmented_frames
            frame_index = cell_frames(segmented_frame_index);
            
            % current_segmentation = cell_segmentations_aligned{frame_index};
            current_segmentation = cell_segmentations_aligned{segmented_frame_index};
            
            % current_translation = [template_parameters.yc, template_parameters.xc, template_parameters.zc] - synapse_centers(frame_index, :);
            % current_translation = [template_parameters.yc, template_parameters.xc, template_parameters.zc] - synapse_centers_smoothed(frame_index, :);
            current_translation = [template_parameters.yc, template_parameters.xc, template_parameters.zc] - synapse_centers_smoothed(segmented_frame_index, :);
            % % Use the centroid instead of the synapse center to prevent issues like with current (r314) version of region_filename = /projects/cellorganizer/xruan/tcell_project/tcell/bhcho_bin_Tcell2/Results/Features_New_Initial_Set3_127_test_snake_regions/CPalpha1/052709 CPalpha1/Run 2/T00045C01Z.mat, per_image_synapse_index = 2:
            % warning('>>>> HACK, have not verified this is correct')
            % current_translation = template_image_centroid - centroids(segmented_frame_index, :);
            % n_x = synapse_normals(frame_index, 1);
            % n_y = synapse_normals(frame_index, 2);
            % n_z = synapse_normals(frame_index, 3);
            % n_x = synapse_normals_smoothed(frame_index, 1);
            % n_y = synapse_normals_smoothed(frame_index, 2);
            % n_z = synapse_normals_smoothed(frame_index, 3);
            n_x = synapse_normals_smoothed(segmented_frame_index, 1);
            n_y = synapse_normals_smoothed(segmented_frame_index, 2);
            n_z = synapse_normals_smoothed(segmented_frame_index, 3);
            [phi, theta, r] = cart2sph(n_x, n_y, n_z);
            theta_degrees = theta * 180 / pi;
            phi_degrees = phi * 180 / pi;

            az = phi_degrees;
            el = theta_degrees;
            current_rotation = [az, el];
            % rotation_center = synapse_centers(frame_index, :)
            % rotation_center = synapse_centers_smoothed(frame_index, :);
            rotation_center = synapse_centers_smoothed(segmented_frame_index, :);
            % % Use the centroid instead of the synapse center
            % rotation_center = centroids(segmented_frame_index, :);
            % point_transform_function = @(given_points)rigidly_transform_points(given_points, maximum_region_size,; current_rotation, current_translation, rotation_center);
            % point_transform_function = @(given_points)rigidly_transform_points(given_points, maximum_region_size, current_rotation, current_translation, rotation_center, true);
            % Now include scaling by synapse radius:
            % current_scale = 1;
            % current_scale = (template_parameters.yr + template_parameters.zr) / 2 / synapse_radii_smoothed(frame_index);
            current_scale = (template_parameters.yr + template_parameters.zr) / 2 / synapse_radii_smoothed(segmented_frame_index);
            point_transform_function_without_scaling = @(given_points)rigidly_transform_points(given_points, maximum_region_size, current_rotation, current_translation, rotation_center, true);
            % point_transform_function = @(given_points)(rigidly_transform_points(given_points, maximum_region_size, current_rotation, current_translation, rotation_center, true) - repmat(rotation_center, [size(given_points, 1), 1])) .* current_scale + repmat(rotation_center, [size(given_points, 1), 1]);
            % point_transform_function = @(given_points)(point_transform_function_without_scaling(given_points) - repmat(rotation_center, [size(given_points, 1), 1])) .* current_scale + repmat(rotation_center, [size(given_points, 1), 1]);
            scaling_transform = translation_matrix3h((rotation_center)) * scale_matrix3h(current_scale .* ones(1, 3)) * translation_matrix3h(-(rotation_center));
            scaling_transform_affine = translation_matrix3h((rotation_center - 1)) * scale_matrix3h(current_scale .* ones(1, 3)) * translation_matrix3h(-(rotation_center - 1));
            % Should create function transform3h to transform 3D points using a matrix for homogeneous coordinates:
            point_transform_function = @(given_points)subsref((scaling_transform * [point_transform_function_without_scaling(given_points), ones(size(given_points, 1), 1)]')', struct('type', '()', 'subs', {{':', 1:3}}));
            
            current_segmentation.vertices = point_transform_function(current_segmentation.vertices);
            current_synapse_point = point_transform_function(current_synapse_point);
            % cell_synapse_points_aligned(frame_index, :) = current_synapse_point;
            cell_synapse_points_aligned(segmented_frame_index, :) = current_synapse_point;
            % segmentation_images{frame_index} = segmentation_rasterization_function(current_segmentation);
            segmentation_images{segmented_frame_index} = segmentation_rasterization_function(current_segmentation);

            % Accumulate transformations with respect to cropped original image:
            % [~, transform, transform_affine] = point_transform_function(zeros(1, 3));
            [~, transform, transform_affine] = point_transform_function_without_scaling(zeros(1, 3));
            % keyboard
            cell_transformation_tracks{current_cell_index, frame_index} = scaling_transform * transform * cell_transformation_tracks{current_cell_index, frame_index};
            cell_transformation_affine_tracks{current_cell_index, frame_index} = scaling_transform_affine * transform_affine * cell_transformation_affine_tracks{current_cell_index, frame_index};
            
            % cell_segmentations_aligned{frame_index} = current_segmentation;
            cell_segmentations_aligned{segmented_frame_index} = current_segmentation;
            % Replace the smoothed variables for degugging visualization below:
            % synapse_normals_smoothed(frame_index, :) = rigidly_transform_points(synapse_normals_smoothed(frame_index, :), maximum_region_size, current_rotation, [0, 0, 0], [0, 0, 0]);
            % synapse_normals_smoothed(frame_index, :) = rigidly_transform_points(synapse_normals_smoothed(frame_index, :), maximum_region_size, current_rotation, [0, 0, 0], [0, 0, 0], true);
            synapse_normals_smoothed(segmented_frame_index, :) = rigidly_transform_points(synapse_normals_smoothed(segmented_frame_index, :), maximum_region_size, current_rotation, [0, 0, 0], [0, 0, 0], true);
            % synapse_normals_smoothed(frame_index, :) = synapse_normals_smoothed(frame_index, :);
            % synapse_centers_smoothed(frame_index, :) = point_transform_function(synapse_centers_smoothed(frame_index, :));
            synapse_centers_smoothed(segmented_frame_index, :) = point_transform_function(synapse_centers_smoothed(segmented_frame_index, :));

            region_filename = [snake_region_location, strrep(synapse_image_path, base_options.image_location, ''), '/', synapse_image_list{frame_index}(1:end-3), 'mat']
            if exist(region_filename, 'file')
              a = load(region_filename)
              
              % cell_synapse_point_tracks{current_cell_index, frame_index}
              
              % keyboard
              per_image_synapse_indices_logical = ~cellfun(@isempty, cell_segmentation_tracks(:, frame_index));
              % current_cell_index
              per_image_synapse_index = find(find(per_image_synapse_indices_logical) == current_cell_index)
              original_raw_region = a.current_result.raw_image{per_image_synapse_index};
              % Use YXZ order for affine:
              % cell_raw_images{current_cell_index, frame_index} = affine(original_raw_region, transform_affine([2, 1, 3, 4], [2, 1, 3, 4]), [], false, [], [], 'crop');
              % cell_raw_images{segmented_frame_index} = affine(original_raw_region, transform_affine([2, 1, 3, 4], [2, 1, 3, 4]), [], false, [], [], 'crop');
              cell_raw_images{segmented_frame_index} = affine(original_raw_region, cell_transformation_affine_tracks{current_cell_index, frame_index}([2, 1, 3, 4], [2, 1, 3, 4]), [], false, [], [], 'crop');
            end
            
            if false
              % imshow(reshape_contrast(segmentation_images{frame_index}, -1))
              imshow(reshape_contrast(segmentation_images{segmented_frame_index}, -1))
              pause
            end
          end

          % % Debugging:
          % if false
            % % Confirm that all synapse normals are [1, 0, 0] to high precision:
            % synapse_normals_smoothed
            % synapse_normals_smoothed_norms = sqrt(sum(synapse_normals_smoothed.^2, 2))

            % clf
            % % set(gcf, 'Position', [100, 100, 0, 0] + [0, 0, number_segmented_frames, 1] * 100)
            % set(gcf, 'Position', [100, 100, 0, 0] + [0, 0, number_segmented_frames, 1 + number_segmented_frames * .2] * 100)
            % hold on
            % max_alpha = .6;
            % % Want total alpha to be max_alpha, max_alpha = a + (a + (a + a * (1 - a)) * (1 - a)) * (1 - a)... = sum(a * (1 - a)^(i - 1), i, 1, n) = a * (1 - a^n) / (1 - a)):
            % % Mathematica: Simplify[Solve[a * (1 - a^n) / (1 - a) == b, a], {0 < a < 1, 0 < b < 1, Element[n, Integers], n > 0}]
            % % Doesn't work...
            % % individual_alpha = max_alpha;
            % individual_alpha = max_alpha / number_segmented_frames;
            % patch_colors = pmkmp(number_segmented_frames, 'IsoL');
            % % all_vertices = zeros(0, 3);
            % all_vertices = cell2mat(cellfun(@(x)x.vertices, cell_segmentations_aligned, 'UniformOutput', false));
            % all_vertices_limits = [min(all_vertices, [], 1); max(all_vertices, [], 1)];
            % all_vertices_ranges = diff(all_vertices_limits, 1, 1);
            
            % synapse_centers_to_plot = synapse_centers_smoothed;
            % synapses_normals_to_plot = synapse_normals_smoothed;
            % synapses_normals_to_plot = synapses_normals_to_plot * mean(std(all_vertices));

            % for frame_index = 1:number_segmented_frames
              % current_segmentation = cell_segmentations_aligned{frame_index};
              % % current_segmentation.vertices(:, 1) = current_segmentation.vertices(:, 1) + all_vertices_ranges(1) * (frame_index - 1);
              % current_segmentation.vertices(:, 2) = current_segmentation.vertices(:, 2) + all_vertices_ranges(2) * (frame_index - 1);
              % synapse_centers_to_plot(frame_index, 2) = synapse_centers_to_plot(frame_index, 2) + all_vertices_ranges(2) * (frame_index - 1);
              % % current_segmentation.vertices(:, 3) = current_segmentation.vertices(:, 3) + all_vertices_ranges(3) * (frame_index - 1);
              % % patch(current_segmentation, 'EdgeColor', [1, .4, .1], 'FaceColor', 'none')
              % patch(current_segmentation, 'EdgeColor', 'none', 'FaceColor', patch_colors(frame_index, :))
              % % patch(current_segmentation, 'EdgeColor', 'none', 'FaceColor', [1, .4, .1], 'FaceAlpha', 1 / number_segmented_frames)
              % % patch(current_segmentation, 'EdgeColor', patch_colors(frame_index, :), 'FaceColor', 'none')
              % % all_vertices = [all_vertices; current_segmentation.vertices];
            % end
            % % view(3)
            % view(atan2(mean(synapse_normals_smoothed(:, 2)), mean(synapse_normals_smoothed(:, 1))) * 180 / pi - 90, 90)
            % set(gca, 'xdir', 'normal', 'ydir', 'normal', 'zdir', 'normal')
            % camlight, lighting phong
            % % sh = plot3(cell_synapse_points_aligned(:, 1), cell_synapse_points_aligned(:, 2), cell_synapse_points_aligned(:, 3), 'bx-');
            % % sh2 = plot3(cell_synapse_points_aligned_smoothed(:, 1), cell_synapse_points_aligned_smoothed(:, 2), cell_synapse_points_aligned_smoothed(:, 3), 'rx-');
            % % ch = plot3(centroids(:, 1), centroids(:, 2), centroids(:, 3), 'yx-');
            % % ch2 = plot3(centroids_smoothed(:, 1), centroids_smoothed(:, 2), centroids_smoothed(:, 3), 'cx-');
            % % legend([sh, sh2, ch, ch2], 'synpts orig', 'synpts smooth', 'centroids orig', 'centroids smooth')
            
            % % quiver_scale = 5;
            % quiver_scale = 0;
            
            % sch2 = plot3(synapse_centers_smoothed(:, 1), synapse_centers_smoothed(:, 2), synapse_centers_smoothed(:, 3), 'rx-');
            % % snh2 = quiver3(synapse_centers_smoothed(:, 1), synapse_centers_smoothed(:, 2), synapse_centers_smoothed(:, 3), synapses_normals_to_plot(:, 1), synapses_normals_to_plot(:, 2), synapses_normals_to_plot(:, 3), quiver_scale, 'b');
            % snh2 = quiver3(synapse_centers_to_plot(:, 1), synapse_centers_to_plot(:, 2), synapse_centers_to_plot(:, 3), synapses_normals_to_plot(:, 1), synapses_normals_to_plot(:, 2), synapses_normals_to_plot(:, 3), quiver_scale, 'b');
            % % legend([sch2, snh2], 'syn ctrs', 'syn nrms')
            % hold off

            % axis equal
            % axis tight

            % title([current_title, ' (synapse points up)'], 'Interpreter', 'none')

            % beep
            % pause
          % end
          
          % Morph these segmentations to the template:
          cropped_size = [template_parameters.imx, template_parameters.imy, template_parameters.imz];
          % Set options for T cell image registration:
          registration_options = struct();
          
          % registration_options.filter_radius = 4;
          registration_options.filter_radius = 8;
          % registration_options.filter_radius = 16;
          % registration_options.filter_radius = 32;
          registration_options.kernel_z_radius = registration_options.filter_radius;
          registration_options.maximum_deformation_per_step = [1, 1, 1] * .5;
          registration_options.window_radius = template_parameters.imx;
          % registration_options.convergence_absolute_error = prod(cropped_size) ./ 2e3;
          % registration_options.convergence_absolute_error = prod(cropped_size) ./ 1e3;
          registration_options.convergence_absolute_error = prod(cropped_size) ./ 5e2;
          registration_options.single_sided = true;

          registration_options.save_intermediates = false;
          registration_options.save_stages = false;
          registration_options.save_intermediate_images = true;
          registration_options.always_save_full_intermediates = true;
          registration_options
          
          % error('Not yet implemented')
          
          cell_filename_prefix = [base_filename, strrep(strrep(fileparts(synapse_file), [base_options.synapse_location, '/'], ''), '/', '_'), num2str(current_cell_index, '_cell%05d_')];
          % for frame_index = 1:number_segmented_frames
          for segmented_frame_index = 1:number_segmented_frames
            frame_index = cell_frames(segmented_frame_index);
            
            frame_deformation_filename = [cell_filename_prefix, num2str(frame_index, 'frame%05d_'), 'morph.mat'];
            [~, frame_deformation_filename_without_directory] = fileparts(frame_deformation_filename);
            [can_start, final_name, final_exists] = chunk_start(base_filename, frame_deformation_filename_without_directory);
            
            if can_start
              % Crop the final aligned segmentation image:
              % cropped_segmentation_image = segmentation_images{frame_index};
              cropped_segmentation_image = segmentation_images{segmented_frame_index};
              cropped_segmentation_image = padarray(cropped_segmentation_image, max(cropped_size - size(cropped_segmentation_image), 0), 'post');
              % whos cropped_segmentation_image, pause
              cropped_segmentation_image = cropped_segmentation_image(1:cropped_size(1), 1:cropped_size(2), 1:cropped_size(3));
              % Set the edges to zero so the registration does not take forever:
              edge_width = 2;
              cropped_segmentation_image(1:edge_width, :, :) = 0;
              cropped_segmentation_image(end - edge_width + 1:end, :, :) = 0;
              cropped_segmentation_image(:, 1:edge_width, :) = 0;
              cropped_segmentation_image(:, end - edge_width + 1:end, :) = 0;
              cropped_segmentation_image(:, :, 1:edge_width) = 0;
              cropped_segmentation_image(:, :, end - edge_width + 1:end) = 0;

              % cropped_raw_image = cell_raw_images{frame_index};
              cropped_raw_image = cell_raw_images{segmented_frame_index};
              % cropped_raw_image = cell_raw_images{current_cell_index, frame_index};
              % cropped_size, size(cropped_raw_image)
              cropped_raw_image = padarray(cropped_raw_image, max(cropped_size - size(cropped_raw_image), 0), 'post');
              % whos cropped_raw_image, pause
              cropped_raw_image = cropped_raw_image(1:cropped_size(1), 1:cropped_size(2), 1:cropped_size(3));

              % Crop the final aligned raw image:
              % a = load(strrep(frame_image_filename, snake_filtered_segmentation_location, snake_region_location))
              % cropped_image = 
              
              source = cropped_segmentation_image;
              % target = get_template(template_parameters);
              target = template_image;
              r = Greedy3D_lambda_pre_compressed(...
                source, target, 1, registration_options)
              
              cropped_segmentation_image_deformed = interp3(cropped_segmentation_image, r.source_deformation{1}.get_data(), r.source_deformation{2}.get_data(), r.source_deformation{3}.get_data());
              cropped_raw_image_deformed = interp3(cropped_raw_image, r.source_deformation{1}.get_data(), r.source_deformation{2}.get_data(), r.source_deformation{3}.get_data());

              % keyboard
              
              save(final_name, 'cropped_segmentation_image', 'cropped_raw_image', 'cropped_segmentation_image_deformed', 'cropped_raw_image_deformed', 'registration_options', 'r')
              
              chunk_finish(base_filename, frame_deformation_filename_without_directory);
            elseif final_exists
              % Turn this on and run as one job to recompute raw image deformation (in case the raw image is not aligned with the segmentation):
              a = load(final_name);
              r = a.r;
              
              % recompute_cropped_raw_image_deformed = false;
              recompute_cropped_raw_image_deformed = true;
              if recompute_cropped_raw_image_deformed
                cropped_raw_image = cell_raw_images{segmented_frame_index};
                cropped_raw_image = padarray(cropped_raw_image, max(cropped_size - size(cropped_raw_image), 0), 'post');
                cropped_raw_image = cropped_raw_image(1:cropped_size(1), 1:cropped_size(2), 1:cropped_size(3));
                cropped_raw_image_deformed = interp3(cropped_raw_image, r.source_deformation{1}.get_data(), r.source_deformation{2}.get_data(), r.source_deformation{3}.get_data());

                a.cropped_raw_image = cropped_raw_image;
                a.cropped_raw_image_deformed = cropped_raw_image_deformed;
                
                save(final_name, '-struct', 'a')
              end
              
              % Create a figure convincing us that the raw image is aligned with the segmentation:
              evaluation_image = reshape_2d([a.cropped_segmentation_image; contrast_stretch(a.cropped_raw_image); contrast_stretch(a.cropped_raw_image) .* a.cropped_segmentation_image; contrast_stretch(a.cropped_raw_image_deformed) .* a.cropped_segmentation_image], -1);
              imwrite(evaluation_image, strrep(final_name, '.mat', '_evaluation.png'))
              if false
                imshow(evaluation_image), pause
              end
              
            end
          end
          

          
            
          % keyboard
          
        % end
      end


      if show_debug_plot
        hold off
        axis equal
        axis tight
        light_handle = lightangle(0, 60); set(light_handle, 'Color', .75 * ones(1, 3));
        light_handle = lightangle(150, -30); set(light_handle, 'Color', .5 * ones(1, 3));
        lighting phong
        colormap(pmkmp(256, 'IsoL'));
        colorbar
        
        % The last frames are not adjacent to each other! This suggests we were expected to connect the dots temporally...
        % synapse_locations_with_frame_indices = cell2mat(arrayfun(@(index)[synapse_locations{index}, ones(size(synapse_locations{index}, 1), 1) * index], 1:length(synapse_locations), 'UniformOutput', false)')

        
        pause
      end

      % keyboard
    end
    
    % % Now need to decide on best one (plot errors relative to landmark_relative_cost first):
    % error('Implementation yet unfinished below this line!')

    error('Implementation yet unfinished below this line!')
    
  end
  
  
  
  
  
  if use_profiling
    p = profile('info');
    profile_filename = [...
      base_filename, ...
      '_profile_results', '.mat']
    profsave(p, profile_filename);
    save([profile_filename, '.mat'], 'p');
  end

  time_so_far = toc(start_time)

  diary('off')
  
end
  
  
