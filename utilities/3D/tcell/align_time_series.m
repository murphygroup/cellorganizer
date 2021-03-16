function [transforms] = align_time_series(options)
  % Rigidly rotate each frame of a time series so the errors between frames and optionally landmarks associated with them are minimized.
  % 
  % To do:
  % Add temporal regularization.
  % 
  % 2013-03-31 tebuck: Created.
  % 2013-04-05 tebuck: Penalizing total intensity reduction to prevent the optimization from moving cells outside the image to reduce total error. This assumes that the pixels have nonnegative values.
  % 10:42 2013-04-10: Added option return_final_errors to get separate cost contributions after optimization finishes.
  
  default_options.number_images = 0;
  default_options.image_function = @(image_index)[];
  default_options.landmark_function = @(image_index)ones(0, 3);
  default_options.image_relative_cost = 1;
  default_options.landmark_relative_cost = 1;
  default_options.absolute_error_image_filename_prefix = '';
  default_options.use_gradient_error_normalization = false;
  % default_options.use_gradient_error_normalization = true;
  default_options.debug = 0;
  % default_options.debug = 1;
  default_options.maximum_function_evaluations = 2000;
  default_options.maximum_iterations = 400;
  % default_options.use_quadratic_cost = true;
  default_options.use_quadratic_cost = false;
  default_options.return_final_errors = false;

  if ~exist('options', 'var')
    options = default_options;
  else
    options = process_options_structure(default_options, options);
  end
  
  options
  
  % Units are voxels:
  % current_translations = zeros(options.number_images, 3);
  current_translations = zeros(options.number_images - 1, 3);
  % Units are degrees:
  % current_rotations = zeros(options.number_images, 2);
  current_rotations = zeros(options.number_images - 1, 2);
  % % Specifies the axis of rotation, where the length of each vector is the number of revolutions:
  % current_rotations = zeros(options.number_images, 3);
  
  function [result_parameter_vector] = parameters_to_parameter_vector(given_translations, given_rotations)
    result_parameter_vector = [given_translations(:); given_rotations(:)];
  end
  
  function [result_translations, result_rotations] = parameter_vector_to_parameters(given_parameter_vector)
    result_translations = reshape(given_parameter_vector(1:numel(current_translations)), size(current_translations));
    result_rotations = reshape(given_parameter_vector(numel(current_translations) + (1:numel(current_rotations))), size(current_rotations));
  end
  
  % function [result_cost] = alignment_cost_function(given_parameter_vector)
  function [varargout] = alignment_cost_function(given_parameter_vector, return_separate_errors)
    % Return the sum of squared errors between consecutive images and landmarks:
    [given_translations, given_rotations] = parameter_vector_to_parameters(given_parameter_vector);
    % result_cost = 0;

    if options.debug > 1
      given_translations, given_rotations
    end
    
    if nargin < 2
      return_separate_errors = false;
    end

    % Do not move the first image:
    first_image = options.image_function(1);
    first_landmarks = options.landmark_function(1);
    previous_image = first_image;
    previous_landmarks = first_landmarks;
    image_result_cost = 0;
    landmark_result_cost = 0;
    
    current_image = first_image;
    for image_index = 2:options.number_images
      if options.image_relative_cost > 0
        current_image = options.image_function(image_index);
        current_image_intensity = sum(current_image(:));
      end
      current_landmarks = options.landmark_function(image_index);
      
      
      % Images are required for centering for now.
      
      if image_index > 1
        % Cumulative transformations, first image is static:
        cumulative_rotation = sum(given_rotations(1:image_index - 1, :), 1);
        cumulative_translation = sum(given_translations(1:image_index - 1, :), 1);
        if options.image_relative_cost > 0
          current_image = rigidly_transform_image(current_image, cumulative_rotation, cumulative_translation);
        end
        if numel(current_landmarks) > 0
          current_landmarks = rigidly_transform_points(current_landmarks, size(first_image), cumulative_rotation, cumulative_translation);
        end
      end
      
      % Compute cost contribution:
      if options.image_relative_cost > 0
        if options.use_gradient_error_normalization
          [~, ~, ~, gradient_magnitude] = normalized_smooth_gradient(first_image);
          error_normalization = sum(gradient_magnitude(:));
        else
          error_normalization = numel(first_image);
        end
        absolute_error = abs(current_image(:) - previous_image(:));
        % Penalizing total intensity reduction the same as misalignment to prevent the optimization from moving cells outside the image to reduce total error:
        absolute_error = [absolute_error; abs(current_image_intensity - sum(current_image(:)))];
        if options.use_quadratic_cost
          image_result_cost = image_result_cost + sum((absolute_error ./ error_normalization).^2);
        else
          image_result_cost = image_result_cost + sum(absolute_error ./ error_normalization);
        end
      end
      
      % Transform current_landmarks by current parameters:
      if numel(current_landmarks) > 0
        % Compute cost contribution:
        error_normalization = numel(first_landmarks) .* mean(std(first_landmarks, 1));
        absolute_error = abs(current_landmarks(:) - previous_landmarks(:));
        if options.use_quadratic_cost
          landmark_result_cost = landmark_result_cost + sum((absolute_error ./ error_normalization).^2, 1);
        else
          landmark_result_cost = landmark_result_cost + sum(absolute_error ./ error_normalization, 1);
          
          % landmark_result_cost_change = sum(absolute_error ./ error_normalization, 1);
          % landmark_result_cost_change
          % if landmark_result_cost_change >= 1e1
            % warning('landmark_result_cost_change >= 1e1!')
            % image_index, absolute_error, error_normalization, landmark_result_cost_change
          % end
        end
      end
      
      if options.debug > 2
        imshow(reshape_contrast([previous_image; current_image])), pause
      end
      
      if options.image_relative_cost > 0
        previous_image = current_image;
      end
      previous_landmarks = current_landmarks;
    end
    
    % result_cost = image_result_cost .* options.image_relative_cost + landmark_result_cost .* options.landmark_relative_cost;
    image_result_cost = image_result_cost .* options.image_relative_cost;
    landmark_result_cost = landmark_result_cost .* options.landmark_relative_cost;
    result_cost = image_result_cost + landmark_result_cost;

    if options.debug > 0
      % image_index, result_cost
      % fprintf('Image %d, result_cost %e\n', image_index, result_cost)
      fprintf('result_cost %e = image_result_cost %e + landmark_result_cost %e\n', result_cost, image_result_cost, landmark_result_cost)
    end
    
    if return_separate_errors
      varargout = {image_result_cost, landmark_result_cost};
    else
      varargout = {result_cost};
    end
  end
  
  
  
  function [should_stop] = debug_plot_misalignment(given_parameter_vector, optimization_values, optimization_state)
    should_stop = false;
    % Return the sum of squared errors between consecutive images and landmarks:
    [given_translations, given_rotations] = parameter_vector_to_parameters(given_parameter_vector);
    % given_translations, given_rotations

    % Do not move the first image:
    first_image = options.image_function(1);
    first_landmarks = options.landmark_function(1);
    previous_image = first_image;
    previous_landmarks = first_landmarks;
    mean_absolute_difference_image = zeros(size(previous_image));
    % Store images as cell entries first:
    images_to_show = cell(1, options.number_images * 2);
    for image_index = 1:options.number_images
      current_image = options.image_function(image_index);
      current_landmarks = options.landmark_function(image_index);
      
      if image_index > 1
        % Cumulative transformations, first image is static:
        cumulative_rotation = sum(given_rotations(1:image_index - 1, :), 1);
        cumulative_translation = sum(given_translations(1:image_index - 1, :), 1);
        current_image = rigidly_transform_image(current_image, cumulative_rotation, cumulative_translation);
        current_landmarks = rigidly_transform_points(current_landmarks, size(current_image), cumulative_rotation, cumulative_translation);
      end
      
      mean_absolute_difference_image = mean_absolute_difference_image + abs(current_image - previous_image);
      
      images_to_show{image_index * 2 - 1} = repmat(current_image, [1, 1, 1, 3]);
      if image_index > 1
        images_to_show{image_index * 2 - 2} = repmat(abs(current_image - previous_image), [1, 1, 1, 3]);
      end
      
      landmark_size = 2;
      rounded_landmarks = round(current_landmarks);
      for landmark_index = 1:size(rounded_landmarks, 1)
        current_rounded_landmark = rounded_landmarks(landmark_index, :);
        if any(current_rounded_landmark < landmark_size) || any(current_rounded_landmark > size(current_image) - landmark_size)
          continue
        end
        images_to_show{image_index * 2 - 1}(current_rounded_landmark(2) + (-1:1), current_rounded_landmark(1) + (-1:1), current_rounded_landmark(3) + (-1:1), 1) = 1;
        images_to_show{image_index * 2 - 1}(current_rounded_landmark(2) + (-1:1), current_rounded_landmark(1) + (-1:1), current_rounded_landmark(3) + (-1:1), 2) = 0;
        images_to_show{image_index * 2 - 1}(current_rounded_landmark(2) + (-1:1), current_rounded_landmark(1) + (-1:1), current_rounded_landmark(3) + (-1:1), 3) = 0;
      end

      previous_image = current_image;
      previous_landmarks = current_landmarks;
    end
    
    % % First frame:
    % images_to_show{1} = repmat(first_image, [1, 1, 1, 3]);

    % Total error:
    images_to_show{end} = repmat(mean_absolute_difference_image, [1, 1, 1, 3]);
    
    % images_to_show'
    % image_to_show = cell2mat(cellfun(@(x)cat(3, reshape_2d(x(:, :, :, 1)), reshape_2d(x(:, :, :, 2)), reshape_2d(x(:, :, :, 3))), images_to_show, 'UniformOutput', false));
    image_to_show = cell2mat(cellfun(@(x)cat(3, reshape_2d(x(:, :, :, 1), size(x, 3)), reshape_2d(x(:, :, :, 2), size(x, 3)), reshape_2d(x(:, :, :, 3), size(x, 3))), images_to_show, 'UniformOutput', false));

      % mean_absolute_difference_image = mean_absolute_difference_image ./ options.number_images;
    % image_to_show = reshape_contrast(mean_absolute_difference_image, 4);
    % image_to_show = reshape_contrast(mean_absolute_difference_image, floor(sqrt(size(previous_image, 3))));
    % image_to_show = reshape_2d(mean_absolute_difference_image, floor(sqrt(size(previous_image, 3))));
    if options.debug > 0
      % imshow(image_to_show)
    end
    if ~isempty(options.absolute_error_image_filename_prefix)
      % imwrite(image_to_show, [options.absolute_error_image_filename_prefix, num2str(optimization_values.iteration, '%05d'), '.png'])
      imwrite(image_to_show, [options.absolute_error_image_filename_prefix, num2str(optimization_values.iteration + 1, '%05d'), '.png'])
    end
    % pause
    drawnow
  end
  
  
  initial_parameter_vector = parameters_to_parameter_vector(current_translations, current_rotations);
  % optimized_parameter_vector = fmincon(@alignment_cost_function, initial_parameter_vector);
  % optimized_parameter_vector = fminunc(@alignment_cost_function, initial_parameter_vector, optimoptions('fminunc', 'TolX', .25, 'TypicalX', initial_parameter_vector .* 0 + 2));
  % optimized_parameter_vector = fminunc(@alignment_cost_function, initial_parameter_vector, struct('TolX', .25, 'TypicalX', initial_parameter_vector .* 0 + 2));
  % Does not fully rotate cubes into each other during test:
  % fminunc_options = struct('TolX', 1e-2, 'TypicalX', (initial_parameter_vector .* 0 + 2) ./ sqrt(eps), 'DiffMaxChange', .5, 'DiffMinChange', 5e-2);
  % fminunc_options = struct('TolX', 1e-3, 'TypicalX', (initial_parameter_vector .* 0 + 2e-1) ./ sqrt(eps), 'DiffMaxChange', .5, 'DiffMinChange', 5e-3);
  % Not much apparent difference here:
  % fminunc_options = struct('TolX', 1e-6, 'TolFun', 1e-6, 'TypicalX', (initial_parameter_vector .* 0 + 2e-1) ./ sqrt(eps), 'DiffMaxChange', .5, 'DiffMinChange', 5e-3);
  % fminunc_options = struct('TolX', 1e-9, 'TolFun', 1e-9, 'TypicalX', (initial_parameter_vector .* 0 + 2e-4) ./ sqrt(eps), 'DiffMaxChange', .5, 'DiffMinChange', 5e-3);
  % fminunc_options = struct('TolX', 1e-9, 'TolFun', 1e-9, 'TypicalX', (initial_parameter_vector .* 0 + 2e-2) ./ sqrt(eps), 'DiffMaxChange', .5, 'DiffMinChange', 5e-3);
  fminunc_options = struct('TolX', 1e-9, 'TolFun', 1e-9, 'TypicalX', [2e0 .* ones(numel(current_translations), 1); 1.5e1 .* ones(numel(current_rotations), 1)] * 1e-4 ./ sqrt(eps), 'DiffMaxChange', 10, 'DiffMinChange', 1e-2);
  % fminunc_options = struct('TolX', 1e-4, 'TolFun', 1e-6, 'TypicalX', (initial_parameter_vector .* 0 + 2e-2) ./ sqrt(eps), 'DiffMaxChange', .5, 'DiffMinChange', 5e-3);
  if options.debug > 1
    % fminunc_options.PlotFcns = {@optimplotfval, @debug_plot_misalignment};
    fminunc_options.PlotFcns = {@debug_plot_misalignment};
  else
    fminunc_options.OutputFcn = {@debug_plot_misalignment};
  end
  
  
  fminunc_options.MaxFunEvals = options.maximum_function_evaluations;
  fminunc_options.MaxIter = options.maximum_iterations;
  fminunc_options.Display = 'iter';
  
  % Write the initial state image so we know what it looks like:
  debug_plot_misalignment(initial_parameter_vector, struct('iteration', -1), struct())
  
  optimized_parameter_vector = fminunc(@alignment_cost_function, initial_parameter_vector, fminunc_options);
  % optimized_parameter_vector = lsqnonlin(@(x)alignment_cost_function(x, false), initial_parameter_vector, [], [], fminunc_options);
  % optimized_parameter_vector = lsqnonlin(@(x)alignment_cost_function(x, true), initial_parameter_vector, [], [], fminunc_options);
  [optimized_translations, optimized_rotations] = parameter_vector_to_parameters(optimized_parameter_vector);
  
  transforms = struct;
  transforms.translations = optimized_translations;
  transforms.rotations = optimized_rotations;
  
  if options.return_final_errors
    [transforms.final_image_error, transforms.final_landmark_error] = alignment_cost_function(optimized_parameter_vector, true);
    if transforms.final_image_error > 0 && options.image_relative_cost > 0
      transforms.final_image_error = transforms.final_image_error / options.image_relative_cost;
    end
    if transforms.final_landmark_error > 0 && options.landmark_relative_cost > 0
      transforms.final_landmark_error = transforms.final_landmark_error / options.landmark_relative_cost;
    end
  end
end


