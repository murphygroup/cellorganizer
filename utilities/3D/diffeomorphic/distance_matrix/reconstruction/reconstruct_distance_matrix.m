function [reconstructed_matrix] = reconstruct_distance_matrix(known_distances, desired_dimensionality, options)
  % Directly solve for all entries in a distance matrix that is partially observed. In known_distances, entries are NaN when unknown. At least desired_dimensionality + 2 columns are assumed to be completely observed.
  %
  % We know of no theoretical guarantees of correctness, accuracy, or precision of this method if known_distances has fewer than (desired_dimensionality + 2) known entries per row/column or is corrupted by noise.
  %
  % Reference: Drineas et al. 2006, "Distance Matrix Reconstruction from Incomplete Distance Information for Sensor Network Localization"
  % 
  % 2013-04-22 tebuck: Copied from svd_reconstruct.m.
  % Dependencies:
  % % From the File Exchange: NIFTI_20121012

  default_options = struct();
  default_options.true_distance_matrix = [];
  % default_options.use_all_known_columns = false;
  default_options.use_all_known_columns = true;
  default_options.use_squared_distances = true;
  default_options.zero_negative_results = true;
  default_options.debug = false;

  if ~exist('options', 'var')
    options = default_options;
  else
    options = process_options_structure(default_options, options);
  end
  
  if options.use_squared_distances
    % Operate on the squared distance matrix:
    known_distances = known_distances.^2;
  end
  
  number_points = size(known_distances, 1);
  
  iteration_index = 1;
  
  reconstructed_matrix = known_distances;
  
  while sum(isnan(known_distances(:))) > 0
    fprintf('Iteration %d\n', iteration_index)
    
    known_value_indices = ~isnan(known_distances);
    % known_value_indices_off_diagonal = known_value_indices & (eye(number_points) == 0);
    
    % Determine reference columns:
    known_column_indices = min(known_value_indices, [], 1) == 1;
    known_columns = known_distances(:, known_column_indices);
    % initial_number_known_columns = sum(known_column_indices);
    % number_known_columns = initial_number_known_columns;
    
    changed = false;

    for column_index = find(~known_column_indices)
      % Entries with which to find the solution for this column:
      known_column_entry_indices = known_value_indices(:, column_index);
      unknown_column_entry_indices = find(~known_column_entry_indices);
      known_column_entry_indices = find(known_column_entry_indices);
      % known_column_entry_indices = known_column_entry_indices(known_column_entry_indices ~= column_index);
      known_column_entries = known_distances(known_column_entry_indices, column_index);
      
      predictors = known_columns(known_column_entry_indices, :);
      predicted = known_column_entries;
      reconstructors = known_columns(unknown_column_entry_indices, :);
      
      if options.use_all_known_columns
        % Determine reference columns for these rows specifically as there may be more information available:
        all_predictor_column_indices = min(known_value_indices(unknown_column_entry_indices, :), [], 1) == 1;
        % number_known_entry_columns = sum(known_entry_column_indices);
        
        % Modify predictors so the reconstruction has more information (predict with knowns, reconstruct the unknowns):
        predictors = known_distances(known_column_entry_indices, all_predictor_column_indices);
        all_predictor_row_indices = min(~isnan(predictors), [], 2) == 1;
        % predictors, all_predictor_row_indices, pause
        predictors = predictors(all_predictor_row_indices, :);
        predicted = known_distances(known_column_entry_indices, column_index);
        predicted = predicted(all_predictor_row_indices, :);
        reconstructors = known_distances(unknown_column_entry_indices, all_predictor_column_indices);
      end

      % predicted, predictors
      % try
        combination_coefficients = pinv(predictors) * predicted;
        % combination_coefficients = predictors \ predicted;
        % Use nonnegative least-squares to ensure we get no negative distances:
        % combination_coefficients = lsqnonneg(predictors, predicted);
        % [combination_coefficients, residual_norm] = lsqnonneg(predictors, predicted);
      % catch err
        % keyboard
        % rethrow(err)
      % end
      
      % keyboard
      % new_values = known_columns(unknown_column_entry_indices, :) * combination_coefficients;
      % new_values = unknown_columns_entries * combination_coefficients;
      % new_values = predictors * combination_coefficients;
      new_values = reconstructors * combination_coefficients;
      % whos new_values
      
      if any(~isfinite(new_values(:)))
        error('Non-finite results!')
      end
      if any(new_values < 0)
        if options.zero_negative_results
          new_values(new_values < 0) = 0;
        end
      
        % beep
        % keyboard
      end
      % whos reconstructed_matrix unknown_column_entry_indices column_index new_values
      reconstructed_matrix(unknown_column_entry_indices, column_index) = new_values;
      reconstructed_matrix(column_index, unknown_column_entry_indices) = reconstructed_matrix(unknown_column_entry_indices, column_index)';
      
      % pause
      % options.debug
      if options.debug
        fprintf('options.debug!\n')
        % pause
        predictors, predicted, combination_coefficients
        % residual_norm
        actual_without = cmdscale(options.true_distance_matrix(known_column_entry_indices, known_column_entry_indices));
        actual_with = cmdscale(options.true_distance_matrix([known_column_entry_indices; unknown_column_entry_indices], [known_column_entry_indices; unknown_column_entry_indices]));
        reconstructed_without = cmdscale(known_distances(known_column_entry_indices, known_column_entry_indices));
        % known_column_entry_indices, unknown_column_entry_indices
        reconstructed_with = cmdscale(reconstructed_matrix([known_column_entry_indices; unknown_column_entry_indices], [known_column_entry_indices; unknown_column_entry_indices]));
        % [~, reconstructed_with_transformed] = procrustes(reconstructed_without, reconstructed_with(1:length(known_column_entry_indices), :));
        % [~, reconstructed_with_transformed] = procrustes(reconstructed_with_known, reconstructed_with(1:length(known_column_entry_indices), :));
        % [~, reconstructed_with_transformed] = procrustes(actual_with, reconstructed_with(1:length(known_column_entry_indices), :));
        [~, reconstructed_with_transformed] = procrustes(actual_with, reconstructed_with);
        clf
        hold on
        % scatter(reconstructed_without(:, 1), reconstructed_without(:,2), 'rx')
        scatter(actual_with(:, 1), actual_with(:, 2), 'rx')
        scatter(reconstructed_with_transformed(:, 1), reconstructed_with_transformed(:, 2), 'bo')
        hold off
        legend({'Actual', 'Reconstructed'})
        pause
        % keyboard
      end
      
    end
    
    if any(known_distances(:) ~= reconstructed_matrix(:))
      changed = true;
    else
      error('Unsuccessful in completing distance matrix')
    end
    
    % reconstructed_matrix(reconstructed_matrix < 0) = nan;
    % reconstructed_matrix = known_distances;
    known_distances = reconstructed_matrix';
    
    % keyboard
    
    iteration_index = iteration_index + 1;
  end
  
  if options.use_squared_distances
    reconstructed_matrix = sqrt(reconstructed_matrix);
  end
  
  
end
