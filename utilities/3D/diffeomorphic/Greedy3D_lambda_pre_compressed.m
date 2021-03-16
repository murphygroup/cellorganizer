function [r_structure] = Greedy3D_lambda_pre_compressed(source, target, lambda_list, options)
  % Greedy3D_lambda_pre_compressed is an LDDMM-based image registration and interpolation method with a wide range of options for control over the image deformation process. Many of these are historical, but possibly still useful, and became necessary at one point or another in order to regulate accuracy, speed, or memory consumption. See Greedy3D_compressed_default_options.m for descriptions of options. Be forewarned that not all combinations of parameters are valid or have been tested.
  % 
  % Arguments:
  % source, 2D or 3D array, the moving image (registration) or first/left/source image (interpolation)
  % target, 2D or 3D array, the fixed image (registration) or second/right/target image (interpolation)
  % lambda_list, 1D array, specifying where between the source and target the interpolated image should be. This can be a ratio, i.e., (distance between interpolated image and target image) / (distance between interpolated image and source image), with 0 being the target, 1 being halfway between them in terms of the LDDMM distance metric, and inf being the source. With options.lambdas_are_distances_from_source = true, each entry can be a distance from the source image. Note that lambda_list can be an array has mostly been used, and tested, as a scalar. Ignored during registration (options.single_sided = true) or distance computation (options.just_compute_distance = true).
  % options, structure, see Greedy3D_compressed_default_options.m for descriptions of options
  %
  % 2011-06-18 tebuck: Copied from Greedy3D_lambda_pre_noneulerian_restarting2.m.
  % 2012-09-20 tebuck: modifying to use WindowedImage for optional in-memory compression.
  % 2013-02-10 tebuck: adding options.lambdas_are_distances_from_source = true for specifying distance from source in lambda_list. This removes an ambiguity when the distances are given as negative values and this value was zero, as it would be treated as a lambda value and so the target was returned instead of the source.
  % 2013-02-17 tebuck: Changing the default options to those that are currently most useful or are being used in CellOrganizer.
  
  %tic
  start_wall_time = tic;
  start_cpu_time = cputime;

  [M, N, P] = size(source);

  % Default option values and options processing (these are now set in Greedy3D_compressed_default_options.m):
  % % known_distance says how far the shapes are known/assumed to be (negative values disable this function):
  % default_options.known_distance = -1;
  % % stop_source_distance is how far the source can move before stopping, and (known_distance - stop_source_distance) is how far the target can move before stopping (negative values disable this function):
  % default_options.stop_source_distance = -1;
  % default_options.single_sided = false; 
  % default_options.filter_radius = 32;
  % default_options.window_radius = 64;
  % default_options.use_compression = false; 
  % default_options.just_compute_distance = false; 
  % default_options.save_intermediates = true; 
  % default_options.lambdas_are_distances_from_source = false; 
  default_options = Greedy3D_compressed_default_options();
  if ~exist('options', 'var')
    options = default_options;
  else
    options = process_options_structure(default_options, options);
  end

  window_size = options.window_radius * [1, 1, 0] + [0, 0, P];
  compression_method = 'none';
  if options.use_compression
    compression_method = 'gzip';
  end
  source_compressed = WindowedImage(source, window_size, compression_method);
  target_compressed = WindowedImage(target, window_size, compression_method);

  clear source target

  % r = Greedy3D_integration_compressed(...
    % source_compressed, target_compressed, ...
    % options.filter_radius, [0, options.max_time], options);
  r = Greedy3D_integration_compressed(source_compressed, target_compressed, options);
  
  
  
  % Number of saved iterations, actually:
  number_iterations = length(r.time_points);
  first_all_t = r.time_points; 
  % Array of distances traveled per iteration:
  distar1 = r.source_distance_array;
  distar2 = [];
  if (~options.single_sided)
    distar2 = r.target_distance_array;
  end
  first_all_distance_array1 = distar1;
  first_all_distance_array2 = distar2;
  % Number of iterations from this first run of Greedy3D_integration_compressed:
  r_structure.total_iterations = r.total_iterations;
  if (options.single_sided)
    % The deformation has been computed and we are done:
    r_structure.interpolated_image = r.source;
    r_structure.total_distance = r.source_distance_array(end); 
    r_structure.source_deformation = r.source_deformation;
    r_structure.total_wall_time = toc(start_wall_time);
    r_structure.total_cpu_time = cputime - start_cpu_time;
    return
  end
  
  %warning('Not finished below this line')
  % Determine the time closest to the one with the correct distance
  % for our lambda:
  first_all_whole_distance = first_all_distance_array1(end) + ...
      first_all_distance_array2(end);
  
  if (length(lambda_list) > 1)
    r_structure.lambda_list = [];
    r_structure.distance_list = [];
    r_structure.interpolated_images = {};
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %for lambda = lambda_list
  for lambda_index = 1:length(lambda_list)
    lambda = lambda_list(lambda_index);
    % fprintf('########################################\n')
    % fprintf('########################################\n')
    % fprintf('lambda_index = %d, lambda = %f:\n\n', lambda_index, lambda)
  
    % The idea is to know that the distance is between two values and
    % also have a guess of the optimal value:
    lambda_time = 0; 
    lambda_time_lower_bound = 0;
    lambda_time_upper_bound = first_all_t(end);
    %lambda_time_lower_bound_index = 1
    lambda_time_lower_bound_index = 0;
    lambda_time_upper_bound_index = length(first_all_t);
    %lambda_side = 0;
    goal_distance = 0;
    goal_left_distance = 0;
    current_distance = 0; 
    %[first_all_t', first_all_distance_array1', first_all_distance_array2']
    from_left_side = true; 
    left_distance = first_all_distance_array1(end);
    right_distance = first_all_distance_array2(end);
    if lambda < 0 || options.lambdas_are_distances_from_source
      goal_distance = lambda;
      if lambda < 0
        warning('LDDMM:lambdaLessThanZero', 'Passing negative distance from source as lambda < 0 is deprecated, use options.lambdas_are_distances_from_source = true instead')
        goal_distance = -lambda;
      end
      % lambda specifies the negative distance to move:
      % 'dist from src 1'
      % goal_distance = -lambda;
      goal_left_distance = goal_distance;
      from_left_side = goal_distance <= left_distance;
    elseif lambda >= first_all_distance_array2(end) / first_all_distance_array1(end)
      % lambda specifies a point within the integration path of the
      % source image:
      % 'lambda from src 1'
      goal_distance = (1/(1+lambda)) * first_all_whole_distance;
      goal_left_distance = goal_distance;
    else
      % lambda specifies a point within the integration path of the
      % target image:
      % 'lambda from trg 1'
      goal_distance = (lambda/(1+lambda)) * first_all_whole_distance;
      goal_left_distance = first_all_whole_distance - goal_distance;
      from_left_side = false;
    end
    % goal_distance
    % goal_left_distance
    % from_left_side

    %wholeDis = distar1(end) + distar2(end);
    wholeDis = first_all_whole_distance;
    final_deformation = []; 
    source_deformation = []; 
    target_deformation = []; 
    ind = 0;
    start_image = [];
    
    img = [];
    if ~options.just_compute_distance
      if options.lambdas_are_distances_from_source && lambda <= 0
        img = source_compressed;
      % If known_distance is not given, use first_all_whole_distance.
      elseif options.lambdas_are_distances_from_source && (lambda >= first_all_whole_distance * (options.known_distance < 0) + options.known_distance * (options.known_distance >= 0))
        img = target_compressed;
      elseif lambda == 0 && ~options.lambdas_are_distances_from_source
        img = target_compressed;
      elseif isinf(lambda)
        img = source_compressed;
      else
        % from_left_side, first_all_t, distar1, distar2, goal_distance, r
        if from_left_side
          [value,ind] = min(abs(distar1 - goal_distance));
          % ind
          start_image = source_compressed; 
          lambda_time = spline(distar1, first_all_t, goal_distance);
        else
          [value,ind] = min(abs(distar2 - goal_distance));
          % ind
          start_image = target_compressed; 
          lambda_time = spline(first_all_whole_distance - distar2, first_all_t, goal_left_distance);
        end
        
        % Linear refinement means linearly interpolate adjacent deformation fields and apply the resulting deformation rather than just the closest one:
        % linear_refinement = ind < r.total_iterations;
        % linear_refinement = ind < length(distar1);
        linear_refinement = (from_left_side && goal_distance < distar1(end)) || (~from_left_side && goal_distance < distar2(end));
        %linear_refinement = false;
        % linear_refinement = true;
        % linear_refinement

        if (linear_refinement)
          % if (lambda_time <= first_all_t(ind))
          if (lambda_time < first_all_t(ind))
            ind = ind - 1;
          end
          if (ind > 0)
            lambda_time_lower_bound = first_all_t(ind);
            lambda_time_lower_bound_index = ind; 
          end
          if (ind < length(first_all_t))
            lambda_time_upper_bound = first_all_t(ind + 1); 
            lambda_time_upper_bound_index = ind + 1;
          end
          % lambda_time
          % lambda_time_lower_bound
          % lambda_time_lower_bound_index
          % lambda_time_upper_bound
          % lambda_time_upper_bound_index
        end


        if options.save_intermediates || (options.known_distance >= 0 && options.stop_source_distance >= 0)
          % r
          if options.save_intermediates
            intermediate_values = load([r.iteration_filenames{ind} '.mat']);
          else
            intermediate_values = r.iteration_steps{ind};
          end
          intermediate_values2 = []; 
          source_deformation = intermediate_values.current_source_deformation;
          target_deformation = intermediate_values.current_target_deformation;

          source_deformation2 = []; 
          target_deformation2 = [];
          final_deformation2 = [];
          if (linear_refinement)
            if options.save_intermediates
              intermediate_values2 = load([r.iteration_filenames{ind + 1} '.mat']);
            else
              % keyboard
              % size(r.iteration_steps), ind
              intermediate_values2 = r.iteration_steps{ind + 1};
            end
            source_deformation2 = intermediate_values2.current_source_deformation; 
            target_deformation2 = intermediate_values2.current_target_deformation; 
          end
          clear intermediate_values intermediate_values2
          if (from_left_side)
            final_deformation = source_deformation; 
            if (linear_refinement)
              final_deformation2 = source_deformation2; 
            end
          else
            final_deformation = target_deformation; 
            if (linear_refinement)
              final_deformation2 = target_deformation2; 
            end
          end

          if (linear_refinement)
            a = (lambda_time - lambda_time_lower_bound) / ...
                (lambda_time_upper_bound - lambda_time_lower_bound)
            %whos
            for dim_ind = 1:3
              final_deformation{dim_ind} = ...
                final_deformation{dim_ind} .* (1 - a) + ...
                final_deformation2{dim_ind} .* a;
            end
          end
        else
          % Recompute the deformation (e.g., when it is faster than reading and writing intermediate files):
          % keyboard
          % r2 = Greedy3D_integration_compressed(...
            % source_compressed, target_compressed, ...
            % options.filter_radius, [0, lambda_time], options);
          options2 = options;
          options2.time_span = [0, lambda_time];
          r2 = Greedy3D_integration_compressed(source_compressed, target_compressed, options2);
          if (from_left_side)
            final_deformation = r2.source_deformation; 
          else
            final_deformation = r2.target_deformation; 
          end
        end

        img = interp3(start_image, final_deformation, '*linear', 0);
      end
    end
      
    
    if (length(lambda_list) > 1)
      r_structure.lambda_list(end + 1)=lambda;
      % Need to compute distance_error at some point:
      %r_structure.distance_list(end + 1)=goal_distance + distance_error;
      r_structure.distance_list(end + 1)=goal_distance;
      r_structure.interpolated_images{end+1} = img;
    else
      r_structure.total_distance = first_all_whole_distance; 
      % Need to compute distance_error at some point:
      %r_structure.deformation_distance = goal_distance + distance_error; 
      r_structure.deformation_distance = goal_distance; 
      if ~options.just_compute_distance
        r_structure.interpolated_image = img; 
        r_structure.source_deformation = source_deformation;
        r_structure.target_deformation = target_deformation;
      end
      r_structure.total_wall_time = toc(start_wall_time);
      r_structure.total_cpu_time = cputime - start_cpu_time;
    end
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end
  
  % Cleanup:
  if ~options.just_compute_distance
    for ind = 1:length(r.iteration_filenames)
      % delete([r.iteration_filenames{ind} '.mat'])
      % Resave much smaller but useful diagnostic files:
      a = load([r.iteration_filenames{ind} '.mat']);
      a.absolute_error = sum(abs(a.current_source - a.current_target));
      a.registration_error = sum((a.current_source - a.current_target).^2);
      % a.maximum_source_deformation_x = max(abs(temp.current_source_deformation{1} - current_source_deformation{1}));
      % a.maximum_source_deformation_y = max(abs(temp.current_source_deformation{2} - current_source_deformation{2}));
      % a.maximum_source_deformation_z = max(abs(temp.current_source_deformation{3} - current_source_deformation{3}));
      % a.maximum_target_deformation_x = max(abs(temp.current_target_deformation{1}) - current_target_deformation{1});
      % a.maximum_target_deformation_y = max(abs(temp.current_target_deformation{2}) - current_target_deformation{2});
      % a.maximum_target_deformation_z = max(abs(temp.current_target_deformation{3}) - current_target_deformation{3});
      a = rmfield(a, {'current_source', 'current_target', 'current_source_deformation', 'current_target_deformation'});
      save([r.iteration_filenames{ind} '.mat'], '-struct', 'a');
    end
  end

  
  if (length(lambda_list) > 1)
    r_structure.total_distance = first_all_whole_distance; 
    r_structure.total_wall_time = toc(start_wall_time);
    r_structure.total_cpu_time = cputime - start_cpu_time;
  end
  
  

  
%finish_wall_time = toc(start_wall_time);

