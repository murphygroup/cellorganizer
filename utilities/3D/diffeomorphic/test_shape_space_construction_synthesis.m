function test_shape_space_construction_synthesis()
  % Trains a shape space model using artificial shapes and then synthesizes novel shapes to demonstrate usage of train_shape_space_model.m and generate_frame_from_shape_space_model.m. Other convenience functions of interest: illustrate_shape_space.m.
  %
  % Dependencies
  % ------------
  % File Exchange: DataHash, isocontour_version2, export_fig
  % Taraz's LDDMM files: 
  % Taraz's shape space files:
  % Rohde's, Tao's, or Wei's shape space files:
  % Taraz's convnfft_fast:
  % 
  % 
  % 2013-02-13 tebuck: Copied from test_generate_frame_from_trained_shape_space_model.m.
  % 2013-06-23 tebuck: Modifying to demonstrate distance matrix completion method. Using a variety of shapes distributed in a three-parameter space but reconstructing a 2D shape space.

  set(0, 'DefaultLineLineWidth', 3, 'DefaultLineMarkerSize', 15)
  set(0, 'DefaultSurfaceLineWidth', 2)
  set(0, 'DefaultAxesFontName', 'Helvetica')
  % Matches my Word font size:
  set(0, 'DefaultAxesFontSize', 11)
  % set(0, 'DefaultAxesFontSize', 12)

  addpath(genpath('DataHash'));
  addpath(genpath('isocontour_version2'));
  addpath(genpath('export_fig'));
  addpath(genpath('pmkmp'));
  % This used to be 'Taraz':
  addpath(genpath('peng_interpolation_code'));

  % I like to use mfilename and parameter identifiers so that tests
  % do not overwrite the results of previous tests:
  filename = mfilename;
  base_filename = [filename, '/'];

  warning('off', 'MATLAB:MKDIR:DirectoryExists');
  mkdir(filename)

  n = getenv('PBS_JOBID');
  n = regexprep(n, '[\r\n]', '');
  date_text = datestr(now(), 'yyyymmddHHMMSSFFF');
  diary([base_filename, 'log', date_text, '_', n, '.txt']);


  use_profiling = false;
%   use_profiling = true;

  if use_profiling
    profile('-memory','on'); profile reset; profile on;
  end
  
  
  start_time = tic;
  start_cputime = cputime;

  try
    
    
    % Either construct/train the shape space directly or from a low-rank approximation to the squared distance matrix (the latter is experimental; we need to find out if this is the least squares way to do it):
    for use_distance_matrix_completion = [false, true]
    
      % current_base_filename = [base_filename];
      % if use_distance_matrix_completion
        % current_base_filename = [base_filename, 'low_rank/'];
      % end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Construct/train the shape space:
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      RandStream.setGlobalStream(RandStream('mt19937ar', 'Seed', 168134));
    
      % Generate some superellipses:
      generate_options = struct();
      generate_options.number_images = 10;
      % generate_options.number_images = 16;
      generate_options.minimum_relative_semidiameter = 1 / 4;
      generate_options.maximum_relative_semidiameter = 2 / 3;
      generate_options.generate_cycle = true;
      generate_options.image_width = 32;
      generate_options
      [shapes, shape_parameters] = generate_superellipse_shape_set(generate_options);
      number_images = length(shapes);
      image_width = size(shapes{1}, 1);

      current_base_filename = [base_filename, sprintf('image_width%05d', image_width)];
      if use_distance_matrix_completion
        current_base_filename = [current_base_filename, '_low_rank'];
      end
      current_base_filename = [current_base_filename, '/'];

      shape_space_options = struct();
      shape_space_options.save_location = [current_base_filename];
      shape_space_options.image_function = @(index)shapes{index};
      shape_space_options.number_images = number_images;
      shape_space_options.voxel_size = [1, 1, 1];
      shape_space_options.shape_aspect_ratio = [1, 1, 1];
      shape_space_options.filter_radius = image_width/4;
      shape_space_options.convergence_absolute_error = (mean([generate_options.minimum_relative_semidiameter, generate_options.maximum_relative_semidiameter]) * image_width * 2 * pi) / 10;
      shape_space_options.maximum_deformation_per_step = ones(1, 3) * .25; 
      if use_distance_matrix_completion
        shape_space_options.desired_shape_space_dimensionality = 2;
        % shape_space_options.desired_shape_space_dimensionality = 5;
      end
      shape_space_options
      
      shape_space = train_shape_space_model(shape_space_options);
      shape_space

      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Plot shapes to visualize the shape space:
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      illustrate_options = struct();
      illustrate_options.image_function = shape_space_options.image_function;
      illustrate_options.number_images = shape_space_options.number_images;
      illustrate_options.voxel_size = shape_space_options.voxel_size;
      illustrate_options.shape_space = shape_space;
      illustrate_options.marker_size_scale = 1;
      illustrate_options.individual_marker_size_influence = 0;
      % illustrate_options.individual_marker_size_influence = .5;
      % illustrate_options.individual_marker_size_influence = .9;
      % illustrate_options.individual_marker_size_influence = 1;
      illustrate_options.shape_levels = [0, 1];
      % illustrate_options.use_outlines = true;
      % illustrate_options.shape_space_inlier_quantile = .5;
      % illustrate_options.shape_space_inlier_quantile = .8;
      % illustrate_options.shape_space_inlier_quantile = .9;
      % illustrate_options.use_third_coordinate_for_color = true;
      illustrate_options.shape_parameters = shape_parameters;
      illustrate_options

      set(gcf, 'Visible', 'off');
      screen_dpi = get(0, 'ScreenPixelsPerInch');
      set(gcf, 'Position', [1, 1, 0, 0] + [0, 0, 6, 4] * screen_dpi);
      set(gcf, 'PaperPositionMode', 'auto');
      use_dark_background = false;
      % use_dark_background = true;
      if use_dark_background
        set(gcf, 'Color', 'k');
        set(gca, 'Color', 'none');
        set(gca, 'XColor', 'w');
        set(gca, 'YColor', 'w');
        set(gca, 'ZColor', 'w');
      else
        set(gcf, 'Color', 'w');
      end
      
      illustrate_shape_space(illustrate_options);
      export_fig([current_base_filename, 'illustration.png'], '-png', '-opengl', '-a1', '-nocrop', ['-r' num2str(300)]); 

      % How well does each coordinate of the shape space correlate with the true shape parameters?
      space_coordinate_vs_parameter_correlations = corr(shape_parameters, shape_space.positions)
      

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Synthesize/generate shapes from the shape space:
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      RandStream.setGlobalStream(RandStream('mt19937ar', 'Seed', 168135));

      % number_images_to_synthesize = 0;
      number_images_to_synthesize = 3;
      % number_images_to_synthesize = 5;
      % number_images_to_synthesize = 500;
      % number_images_to_synthesize = 5000;
      
      if use_distance_matrix_completion
        number_images_to_synthesize = 0;
        warning('Not synthesizing images for the shape space computed from the reconstructed distance matrix due to poor shape space performance.')
      end
      
      random_positions = zeros(number_images_to_synthesize, size(shape_space.positions, 2));
      % combine_many_positions = false;
      combine_many_positions = true;
      % Randomly combine few points (for these shapes, visually nearly uniform):
      number_positions_to_combine = 3;
      % number_positions_to_combine = 5;
      % Randomly combine all points (for these shapes, strongly biased towards the center):
      % number_positions_to_combine = number_images;
      for synthesis_index = 1:number_images_to_synthesize
        % Generate a random point inside the convex hull of training shapes's poitions in the shape space:
        if combine_many_positions
          % Randomly combine number_positions_to_combine points:
          random_weights = zeros(number_images, 1);
          random_vertices = randperm(number_images);
          random_vertices = random_vertices(1:number_positions_to_combine);
          random_weights(random_vertices) = rand(number_positions_to_combine, 1);
        else
          % Randomly combine only two points (for these shapes, weakly biased towards the edges with a uniform sample):
          random_weights = zeros(number_images, 1);
          random_vertices = [randi(number_images, 1), randi(number_images - 1, 1)];
          random_vertices(2) = random_vertices(2) + (random_vertices(2) >= random_vertices(1));
          random_weights(random_vertices(1)) = rand;
          % random_weights(random_vertices(1)) = betarnd(2, 2);
          random_weights(random_vertices(2)) = 1 - random_weights(random_vertices(1));
        end
        % Normalize weights:
        random_weights = random_weights ./ sum(random_weights);
        random_position = random_weights' * shape_space.positions;
        random_positions(synthesis_index, :) = random_position;
      end
      
      hold on
      scatter(random_positions(:, 1), random_positions(:, 2), get(0, 'DefaultAxesFontSize')^2, 'kx')
      export_fig([current_base_filename, 'illustration_with_positions_to_synthesize.png'], '-png', '-opengl', '-a1', '-nocrop', ['-r' num2str(300)]); 
      hold off
        
      for synthesis_index = 1:number_images_to_synthesize
        image_filename = [current_base_filename, sprintf('synthesized_shape%05d.png', synthesis_index)];
        if exist(image_filename, 'file')
          continue
        end
        random_position = random_positions(synthesis_index, :);
        % Generate the shape associated with that point:
        generate_options = shape_space_options;
        generate_options.positions = random_position;
        frame_images = generate_frame_from_shape_space_model(generate_options);
        imwrite(contrast_stretch(frame_images{1}.interpolated_image.get_data()), image_filename)
      end
      
    end
    
    
  catch err
    getReport(err, 'extended')
  end
  
  total_time = toc(start_time)
  total_cputime = cputime - start_cputime

  if use_profiling
    p = profile('info');
    profile_filename = [base_filename, 'profile_results', '.mat'];
    profsave(p, profile_filename);
    save([profile_filename, '.mat'], 'p');
  end
  
  diary('off')

end

