function [positions, mass_matrix, embedding_info] = embed_partial_distance_matrix(distances, embedding_info)
% Provides a consistent interface to compute a shape space from a partially observed distance matrix using one of several methods. The observations must have a certain structure for certain methods, and those should be documented here...
%
% distances is a sparse, square, symmetric, nonnegative matrix with positive off-diagonal values except where the distance has not yet been observed (i.e., there should not be identical objects, and any observed should never round to zero).
%
% Methods nystrom-euclidean, nystrom-dissimilarity-clip, nystrom-dissimilarity-orthogonalize[, others?]: Must observe the first options.desired_shape_space_dimensionality + 2 columns.
%   
% 2013-12-05 tebuck: Copied from complete_distance_matrix.m.
% 2013-12-14 tebuck: distances is now assumed to be sparse (not yet fully implemented!).
% 2013-12-15 tebuck: eventually this should allow parallel computations/always work the same way with a given random seed or number stream...
% 2014-03-15 tebuck: Adding method sdm-nystrom-kernel-orthogonalize.
% 2014-03-16 tebuck: Adding method approximate-svd.


    default_options = struct();
    default_options.desired_shape_space_dimensionality = inf;
    % % If this function is given and the method correctly selected, new distances can be computed upon request (this function should take care of things like using chunk_start, and embed_partial_distance_matrix should use its own random number generator instance instead of the default one):
    % % default_options.distance_function = @(given_index)[];
    % default_options.distance_function = [];
    default_options.method = 'direct';
    % Some randomized methods can run multiple times and select the best embedding:
    default_options.number_random_method_restarts = 1;

    % Clip negative eigenvalues and such so the results have a positive definite distance matrix regardless of the input's definiteness:
    default_options.force_positive_definiteness = true;
    % % Use the Schur decomposition to find the eigendecomposition instead of eig:
    % default_options.use_schur_decomposition = false;
    % % default_options.use_schur_decomposition = true;

    if ~exist('embedding_info', 'var')
    embedding_info = default_options; 
    else
    embedding_info = ml_initparam(embedding_info, default_options);
    end

    embedding_info = ml_initparam(embedding_info, struct( ...
            'explained_variance_threshold', 0.90, ...
            'weight_factor', 0, ...
            'replicates', 1, ...
            'statset', statset('Display', 'iter','MaxIter', 3000, 'TolFun', 1e-6), ...
            'Start', 'random'));

    %remove all the points that we don't have anyinformation about
    distances_incomplete_temp = distances;

    distances_incomplete_temp(logical(eye(size(distances_incomplete_temp)))) = nan;
    keepinds_master = ~all(isnan(distances_incomplete_temp),1);
    distances = distances(keepinds_master, keepinds_master);



    % Copy some options to local variables:
    desired_shape_space_dimensionality = embedding_info.desired_shape_space_dimensionality;

    mass_matrix = [];
    
    number_data = size(distances, 1);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Common code:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    given_entries = distances > 0;
    given_entries_with_diagonals = given_entries | logical(speye(number_data));
    given_entries_proportion = full(sum(sum(given_entries))) ./ numel(given_entries);


    % Find landmarks using already observed entries:

    use_fixed_dimensionality_landmarks = false;
    % use_fixed_dimensionality_landmarks = true;
    if use_fixed_dimensionality_landmarks
        % The number of landmarks is chosen based on a Euclidean matrix's rank if its points have dimensionality desired_shape_space_dimensionality (see Drineas et al. 2006, "Distance Matrix Reconstruction from Incomplete Distance Information for Sensor Network Localization"):
        % warning('Should just check for complete columns or submatrices?')
        number_landmarks = desired_shape_space_dimensionality + 2;
        landmark_indices = [true(1, number_landmarks), false(1, number_data - number_landmarks)];
    else
        % Number of observations not on the diagonal (assumes no duplicate points):
        % column_number_observations = sum(~isnan(distances) & distances ~= 0);
        column_number_observations = sum(given_entries);
        landmark_indices = column_number_observations == max(column_number_observations);
        number_landmarks = sum(landmark_indices);
    end
    % keyboard

    if number_landmarks == 0
    error('No landmarks observed! No methods currently complete randomly sampled distance matrices.')
        % for novel_index = randi(number_data)
          % for second_index = [1:novel_index - 1, novel_index + 1:number_data]
            % if distances(novel_index, second_index) ~= 0
              % continue
            % end
            % distances(novel_index, second_index) = options.distance_function(novel_index, second_index);
            % % Prevent zero distances:
            % if distances(novel_index, second_index) == 0
              % distances(novel_index, second_index) = eps;
              % % distances(novel_index, second_index) = eps(mean(mean(distances)));
            % end
            % distances(second_index, novel_index) = distances(novel_index, second_index);
          % end

        % end
    end

    % if all(column_number_observations) == number_data - 1
    % end

  
  function [result_eigenvectors, result_eigenvalues] = spectral_decomposition(given_matrix, given_sort_by_absolute_eigenvalues)
    if ~exist('given_sort_by_absolute_eigenvalues', 'var') || isempty(given_sort_by_absolute_eigenvalues)
      given_sort_by_absolute_eigenvalues = false;
    end
    % if options.use_schur_decomposition
      % [result_eigenvectors, result_eigenvalues] = schur(given_matrix);
    % else
      [result_eigenvectors, result_eigenvalues] = eig(given_matrix);
    % end
    % Because schur sometimes returns non-diagonal second output:
    result_eigenvalues = diag(result_eigenvalues);
    if given_sort_by_absolute_eigenvalues
      % Sort by decreasing absolute eigenvalue:
      [~, result_eigenvalue_order] = sort(abs(result_eigenvalues), 1, 'descend');
    else
      % Sort by decreasing eigenvalue:
      [~, result_eigenvalue_order] = sort(result_eigenvalues, 1, 'descend');
    end
    result_eigenvalues = result_eigenvalues(result_eigenvalue_order);
    result_eigenvectors = result_eigenvectors(:, result_eigenvalue_order);
    % % Because this can be sparse:
    % result_eigenvalues = spdiags(result_eigenvalues, 0, number_landmarks, number_landmarks);
    % No sparsity for now for convenience:
    result_eigenvalues = diag(result_eigenvalues);
  end
  
  
  function [result_square_root] = general_square_root(given_matrix)
    [given_matrix_eigenvectors, given_matrix_eigenvalues] = spectral_decomposition(given_matrix);
    % result_square_root = given_matrix_eigenvectors * sqrtm(given_matrix_eigenvalues) * pinv(given_matrix_eigenvectors);
    given_matrix_eigenvalues_square_root = diag(sqrt(diag(given_matrix_eigenvalues)));
    result_square_root = given_matrix_eigenvectors * given_matrix_eigenvalues_square_root * pinv(given_matrix_eigenvectors);
  end
  
  
  function [result_eigenvectors, result_eigenvalues] = filter_spectral_decomposition(given_eigenvectors, given_eigenvalues, given_sort_by_absolute_eigenvalues)
    if ~exist('given_sort_by_absolute_eigenvalues', 'var') || isempty(given_sort_by_absolute_eigenvalues)
      given_sort_by_absolute_eigenvalues = false;
    end
    result_eigenvalues = given_eigenvalues;
    result_eigenvectors = given_eigenvectors;
    % Because schur sometimes returns non-diagonal second output:
    result_eigenvalues = diag(result_eigenvalues);
    if given_sort_by_absolute_eigenvalues
      % Sort by decreasing absolute eigenvalue:
      [~, result_eigenvalue_order] = sort(abs(result_eigenvalues), 1, 'descend');
    else
      % Sort by decreasing eigenvalue:
      [~, result_eigenvalue_order] = sort(result_eigenvalues, 1, 'descend');
    end
    result_eigenvalues = result_eigenvalues(result_eigenvalue_order);
    result_eigenvectors = result_eigenvectors(:, result_eigenvalue_order);
    % % Because this can be sparse:
    % result_eigenvalues = spdiags(result_eigenvalues, 0, number_landmarks, number_landmarks);
    % No sparsity for now for convenience:
    result_eigenvalues = diag(result_eigenvalues);

    unfiltered_result_eigenvalues = diag(result_eigenvalues).';
    embedding_info.unfiltered_result_eigenvalues = unfiltered_result_eigenvalues;
    
    good_eigenvalues = abs(unfiltered_result_eigenvalues) > eps(max(abs(unfiltered_result_eigenvalues))) * 1e2;
    
    if embedding_info.force_positive_definiteness
      unfiltered_eigenvalues_negative_proportion = abs(sum(unfiltered_result_eigenvalues(unfiltered_result_eigenvalues < 0))) / sum(abs(unfiltered_result_eigenvalues));
      
      good_eigenvalues = good_eigenvalues & (unfiltered_result_eigenvalues > 0);
      
      embedding_info.unfiltered_eigenvalues_negative_proportion = unfiltered_eigenvalues_negative_proportion;
    end

    % Respect desired_shape_space_dimensionality:
    % good_eigenvalues(desired_shape_space_dimensionality + 1:end) = false;
    good_eigenvalues(cumsum(good_eigenvalues) > desired_shape_space_dimensionality) = false;
    
    % Apply eigenvalue filtering:
    result_eigenvalues = diag(unfiltered_result_eigenvalues .* good_eigenvalues);
    
    eigenvalue_reordering = [find(good_eigenvalues), find(~good_eigenvalues)];
    result_eigenvalues = result_eigenvalues(eigenvalue_reordering, eigenvalue_reordering);
    result_eigenvectors = result_eigenvectors(:, eigenvalue_reordering);
    good_eigenvalues = good_eigenvalues(eigenvalue_reordering);
    result_eigenvalues = result_eigenvalues(good_eigenvalues, good_eigenvalues);
    result_eigenvectors = result_eigenvectors(:, good_eigenvalues);

  end
  
  
  % Initialize Mishra et al. 2011 method:
  mishra_parameters = struct();
  % Wasn't even doing anything, no such option rmax!
  % mishra_parameters.rmax = desired_shape_space_dimensionality + 2;
  % mishra_parameters.pmax = desired_shape_space_dimensionality;
  mishra_parameters.pmax = desired_shape_space_dimensionality + 2;
  mishra_parameters.verb = false;
  % From the demo:
  mishra_parameters.tol = 1e-3;
  % % Default:
  % mishra_parameters.tol = 1e-5;
  % % Want relative convergence, guess this should be low:
  % mishra_parameters.tol = 1e-7;
  % From the demo:
  mishra_parameters.vtol = 1e-3;
  mishra_method_function = @tr_dist_completion;
  % Produces bad results with current options 21:05 2013-12-29...
  % mishra_method_function = @gd_dist_completion;
  % mishra_method_function
  
  
  should_stress_use_glimmer_definition = false;
  % should_stress_use_glimmer_definition = true;
  
  function [result_stress] = partial_stress_function(given_positions, given_subsampling_rate)
    if ~exist('given_subsampling_rate', 'var')
      given_subsampling_rate = 1;
    end
    given_entries_to_check = given_entries;
    given_entries_to_keep = rand(size(given_entries_to_check(given_entries)));
    given_entries_to_keep = given_entries_to_keep <= quantile(given_entries_to_keep(:), given_subsampling_rate);
    given_entries_to_check(given_entries) = given_entries_to_check(given_entries) & given_entries_to_keep;
    [given_entries_to_check_rows, given_entries_to_check_columns] = find(given_entries_to_check);
    
    result_stress = 0;
    for entry_index = 1:length(given_entries_to_check_rows)
      row_index = given_entries_to_check_rows(entry_index);
      column_index = given_entries_to_check_columns(entry_index);
      % result_stress = result_stress + (distances(row_index, column_index).^2 - sum(diff(given_positions([row_index, column_index], :), 1), 2)).^2;
      if should_stress_use_glimmer_definition
        % Glimmer definition:
        result_stress = result_stress + (distances(row_index, column_index) - sqrt(sum(diff(given_positions([row_index, column_index], :), 1).^2, 2))).^2;
      else
        result_stress = result_stress + (distances(row_index, column_index).^2 - sum(diff(given_positions([row_index, column_index], :), 1).^2, 2)).^2;
      end
    end
    
    if should_stress_use_glimmer_definition
      % Glimmer definition:
      result_stress = result_stress ./ sum(distances(given_entries_to_check).^2);
    else
      result_stress = result_stress ./ sum(distances(given_entries_to_check).^4);
    end
  end
  
  function [result_stress_gradient] = partial_stress_gradient_function(given_positions, given_subsampling_rate)
    if ~exist('given_subsampling_rate', 'var')
      given_subsampling_rate = 1;
    end
    given_entries_to_check = given_entries;
    given_entries_to_keep = rand(size(given_entries_to_check(given_entries)));
    given_entries_to_keep = given_entries_to_keep <= quantile(given_entries_to_keep(:), given_subsampling_rate);
    given_entries_to_check(given_entries) = given_entries_to_check(given_entries) & given_entries_to_keep;
    [given_entries_to_check_rows, given_entries_to_check_columns] = find(given_entries_to_check);
    
    result_stress_gradient = zeros(size(given_positions));
    for entry_index = 1:length(given_entries_to_check_rows)
      row_index = given_entries_to_check_rows(entry_index);
      column_index = given_entries_to_check_columns(entry_index);
      for coordinate_index = 1:size(given_positions, 2)
        if should_stress_use_glimmer_definition
          % Glimmer definition:
          error('Not yet implemented!')
          % result_stress_gradient = result_stress_gradient - 4 * (diff(given_positions([row_index, column_index], ))) * (distances(row_index, column_index) - sqrt(sum(diff(given_positions([row_index, column_index], :), 1).^2, 2))).^2;
        else
          % error
          % -4 * (coord diff) * (squared dist diff):
          current_row_term = -4 * (given_positions(row_index, coordinate_index) - given_positions(column_index, coordinate_index)) * (distances(row_index, column_index).^2 - sum(diff(given_positions([row_index, column_index], :), 1).^2, 2));
          result_stress_gradient(row_index, coordinate_index) = result_stress_gradient(row_index, coordinate_index) + current_row_term;
          % result_stress_gradient(column_index, coordinate_index) = result_stress_gradient(column_index, coordinate_index) - current_row_term;
        end
      end
    end
    
    if should_stress_use_glimmer_definition
      % Glimmer definition:
      % error
      result_stress_gradient = result_stress_gradient ./ sum(distances(given_entries_to_check).^2);
    else
      % error
      result_stress_gradient = result_stress_gradient ./ sum(distances(given_entries_to_check).^4);
    end
  end
  
  function [result_stress, result_stress_gradient] = partial_stress_with_gradient_function(given_positions, given_subsampling_rate)
    if ~exist('given_subsampling_rate', 'var')
      given_subsampling_rate = 1;
    end
    result_stress = partial_stress_function(given_positions, given_subsampling_rate);
    result_stress_gradient = partial_stress_gradient_function(given_positions, given_subsampling_rate);
    % Because apparently GradObj output doesn't respect matrix shape even though the inputs do:
    result_stress_gradient = result_stress_gradient(:);
  end
  
  % % Mass matrix is the identity by default:
  
  % mass_matrix = eye(number_landmarks);
  
  
  % embedding_info empty by default:
  
%   embedding_info = struct();
  embedding_info.number_landmarks = number_landmarks;
  embedding_info.landmark_indices = landmark_indices;
  
  
  embedding_info.total_wall_time = tic;
  embedding_info.total_cpu_time = cputime;
  
  
  % Reconstruct entries with nans:
  switch embedding_info.method
    case 'direct'
    
    
      % error('Not yet implemented!')
      % if options.force_positive_definiteness
        % error('Not yet implemented!')
      % end
      % warning('Should respect options.desired_shape_space_dimensionality!')
      
      positions = cmdscale(complete_distance_matrix(distances, embedding_info));
      positions = positions(:, 1:min(end, desired_shape_space_dimensionality));
      mass_matrix = eye(size(positions, 2));
      
      
    case 'mdscale'
      
      
      % keyboard
      [positions, stress] = mdscale(distances, desired_shape_space_dimensionality, 'Criterion', 'strain', 'Start', 'random', 'Weights', given_entries);
      embedding_info.stress = stress;
      positions = positions(:, 1:min(end, desired_shape_space_dimensionality));
      mass_matrix = eye(size(positions, 2));
      
      
    case 'mdscale-squared'
      
      
      [positions, stress] = mdscale(distances.^2, desired_shape_space_dimensionality, 'Criterion', 'strain', 'Start', 'random', 'Weights', given_entries);
      embedding_info.stress = stress;
      positions = positions(:, 1:min(end, desired_shape_space_dimensionality));
      mass_matrix = eye(size(positions, 2));
      
      
    case 'stress-minimization'
      
      
      % 2013-12-30 tebuck: Wanted to implement my own stress minimization using missing values, but back to Brand 2002 strategy!
      % error('Not yet implemented!')
      
      % Could make this deal with Klein spaces, too...
      
      potential_position_sets = cell(embedding_info.number_random_method_restarts, 1);
      potential_position_set_stresses = inf(embedding_info.number_random_method_restarts, 1);
      for restart_index = 1:embedding_info.number_random_method_restarts
        positions = randn(size(distances, 1), desired_shape_space_dimensionality);
        % [positions, stress] = fminunc(@partial_stress_function, positions);
        % [positions, stress] = fminunc(@partial_stress_function, positions, optimset('Display', 'iter', 'MaxFunEvals', numel(positions) * 100 * 5, 'TolFun', 1e-6));
        % [positions, stress] = fminunc(@partial_stress_function, positions, optimset('Display', 'iter', 'MaxFunEvals', numel(positions) * 100 * 5, 'TolFun', 1e-3));
        % [positions, stress] = fminunc(@partial_stress_function, positions, optimset('Display', 'iter', 'MaxFunEvals', numel(positions) * 100 * 5, 'TolFun', 1e-3 * given_entries_proportion));
        [positions, stress] = fminunc(@partial_stress_with_gradient_function, positions, optimset('Display', 'iter', 'MaxFunEvals', numel(positions) * 100 * 5, 'TolFun', 1e-3 * given_entries_proportion, 'GradObj', 'on'));
        % positions = fminsearch(@partial_stress_function, positions, optimset('Display', 'iter'));
        potential_position_sets{restart_index} = positions;
        potential_position_set_stresses(restart_index) = stress;
      end
      [~, minimum_stress_index] = min(potential_position_set_stresses);
      positions = potential_position_sets{minimum_stress_index};
      mass_matrix = eye(size(positions, 2));
      
      embedding_info.potential_position_sets = potential_position_sets;
      embedding_info.potential_position_set_stresses = potential_position_set_stresses;
      embedding_info.minimum_stress_index = minimum_stress_index;
      
      
    % case 'stress-minimization-nmds-init'
      
      
      % % 2013-12-30 tebuck: Wanted to implement my own stress minimization using missing values, but back to Brand 2002 strategy!
      % % error('Not yet implemented!')
      
      % % Could make this deal with Klein spaces, too...
      
      % nmds_positions = ;
      % [positions, stress] = fminunc(@partial_stress_with_gradient_function, nmds_positions, optimset('Display', 'iter', 'MaxFunEvals', numel(nmds_positions) * 100 * 5, 'TolFun', 1e-3 * given_entries_proportion, 'GradObj', 'on'));
      % positions = potential_position_sets{minimum_stress_index};
      % mass_matrix = eye(size(positions, 2));
      
      % embedding_info.nmds_positions = nmds_positions;
      % embedding_info.stress = stress;
      
      
    case 'landmark-kernel-svd'
      
      
      % error('Not yet implemented!')
      % error('Not right!')
      % Based on Sakai and Imiya 2009, "Fast Spectral Clustering with Random Projection and Sampling," algorithm 3, and own reasoning, slightly different than method similarity-principal-components below.

      % Obtain A and B (variable names from Platt 2005):
      landmark_distances = full(distances(landmark_indices, landmark_indices));
      remainder_distances = full(distances(~landmark_indices, landmark_indices));
      
      % Center A and B as in Landmark MDS:
      % A, Platt 2005, equation (13):
      landmark_centering_matrix = eye(number_landmarks) - 1 / number_landmarks;
      % This becomes landmark_centering_matrix * [E_ij^2 - 1/m * sum_q E_iq^2] = [(E_ij^2 - 1/m * sum_q E_iq^2) - sum_p (E_pj^2 - 1/m * sum_q E_pq^2)]?
      landmark_kernel = -1 / 2 * landmark_centering_matrix * (landmark_distances.^2) * landmark_centering_matrix;
      % B, Platt 2005, equation (15) (could change this to (14) at some point):
      % remainder_kernel = -1 / 2 * (remainder_distances.^2 - repmat(mean(remainder_distances.^2, 2), [1, number_landmarks]));
      remainder_kernel = -1 / 2 * (remainder_distances.^2 - repmat(mean(landmark_distances.^2, 1), [number_data - number_landmarks, 1]));
      
      % columns_kernel = [landmark_kernel; remainder_kernel];
      columns_kernel = [landmark_kernel; remainder_kernel];
      columns_kernel = zeros(number_data, size(landmark_kernel, 2));
      columns_kernel(landmark_indices, :) = landmark_kernel;
      columns_kernel(~landmark_indices, :) = remainder_kernel;

      [columns_kernel_u, columns_kernel_s, columns_kernel_v] = svd(columns_kernel, 'econ');
      
      positions = columns_kernel_u * sqrt(columns_kernel_s);
      % positions = columns_kernel_u;
      % positions = columns_kernel * columns_kernel_v * sqrt(columns_kernel_s);
      % positions = columns_kernel * columns_kernel_v * pinv(columns_kernel_s);
      % keyboard

      positions = positions(:, 1:min(end, desired_shape_space_dimensionality));
      mass_matrix = eye(size(positions, 2));
      
      
    case 'landmark-distance-svd'
      
      
      % error('Not yet implemented!')
      % error('Not right!')
      % Based on Sakai and Imiya 2009, "Fast Spectral Clustering with Random Projection and Sampling," algorithm 3, and own reasoning, slightly different than method similarity-principal-components below.

      % Obtain A and B (variable names from Platt 2005):
      landmark_distances = full(distances(landmark_indices, landmark_indices));
      remainder_distances = full(distances(~landmark_indices, landmark_indices));
      
      % columns_distances = [landmark_distances; remainder_distances];
      columns_distances = zeros(number_data, size(landmark_distances, 2));
      columns_distances(landmark_indices, :) = landmark_distances;
      columns_distances(~landmark_indices, :) = remainder_distances;
      [columns_u, columns_s, columns_v] = svd(columns_distances, 'econ');
      
      positions = columns_u * sqrt(columns_s);
      % positions = columns_u;
      % positions = columns * columns_v * sqrt(columns_s);
      % positions = columns * columns_v * pinv(columns_s);
      
      if false
        % Debug info:
        set(gcf, 'visible', 'on')
        clf, hold on, scatter(positions(:, 1), positions(:, 2), 'bx'), legend('Censored recovered aligned'), axis image
        keyboard
      end

      positions = positions(:, 1:min(end, desired_shape_space_dimensionality));
      mass_matrix = eye(size(positions, 2));
      
      
    case 'landmark-distance-squared-svd'
      
      
      % error('Not yet implemented!')
      % Based on Sakai and Imiya 2009, "Fast Spectral Clustering with Random Projection and Sampling," algorithm 3, and own reasoning, slightly different than method similarity-principal-components below.

      % Obtain A and B (variable names from Platt 2005):
      landmark_distances = full(distances(landmark_indices, landmark_indices));
      remainder_distances = full(distances(~landmark_indices, landmark_indices));
      
      % columns_squared_distances = [landmark_distances; remainder_distances].^2;
      columns_squared_distances = zeros(number_data, size(landmark_distances, 2));
      columns_squared_distances(landmark_indices, :) = landmark_distances;
      columns_squared_distances(~landmark_indices, :) = remainder_distances;
      columns_squared_distances = columns_squared_distances.^2;
      [columns_u, columns_s, columns_v] = svd(columns_squared_distances, 'econ');
      
      positions = columns_u * sqrt(columns_s);
      % positions = columns_u;
      % positions = columns * columns_v * sqrt(columns_s);
      % positions = columns * columns_v * pinv(columns_s);
      
      if false
        % Debug info:
        set(gcf, 'visible', 'on')
        clf, hold on, scatter(positions(:, 1), positions(:, 2), 'bx'), legend('Censored recovered aligned'), axis image
        keyboard
      end

      positions = positions(:, 1:min(end, desired_shape_space_dimensionality));
      mass_matrix = eye(size(positions, 2));
      
      
    case 'incremental-svd'
      
      
      % 2013-12-30 tebuck: Wanted to implement my own stress minimization using missing values, but back to Brand 2002 strategy!
      % error('Not yet implemented!')
      
      use_my_implementation = false;
      % use_my_implementation = true;
      
      % This assumes observations exist for every row and column and only centers observed distances:
      distances_squared_centered = distances.^2;
      row_means = sum(distances_squared_centered, 2) ./ sum(given_entries_with_diagonals, 2);
      column_means = sum(distances_squared_centered, 1) ./ sum(given_entries_with_diagonals, 1);
      overall_mean = mean(distances_squared_centered(given_entries_with_diagonals));
      [rows, columns] = find(given_entries_with_diagonals);
      for entry_index = 1:length(rows)
        row = rows(entry_index);
        column = columns(entry_index);
        % Platt 2005 (1):
        distances_squared_centered(row, column) = -1/2 * (distances_squared_centered(row, column) - row_means(row) - column_means(column) + overall_mean);
      end
      % distances_squared_centered(given_entries_with_diagonals) = distances_squared_centered(given_entries_with_diagonals) ./ ;
      
      % [u, s, v] = masked_svd(distances_squared_centered, given_entries_with_diagonals);
      % sparse_eigenvectors = u;
      % sparse_eigenvalues = diag(diag(s) .* sign(sum(u .* v, 1)));
      
        
      % % Not really necessary, masked_svd will usually return a fairly high-rank decomposition and can't currently ask it for a higher rank one...
      
      current_dimensionality = 0;
      number_additional_dimensions = 0;
      
      while current_dimensionality < desired_shape_space_dimensionality && (desired_shape_space_dimensionality + 2 + number_additional_dimensions <= number_data)
        
        
        if use_my_implementation
          
          % [u, s, v] = masked_svd(distances_squared_centered, given_entries_with_diagonals);
          [u, s, v] = masked_svd(distances_squared_centered, given_entries_with_diagonals, struct('rank_limit', desired_shape_space_dimensionality + 2 + number_additional_dimensions));
          
        else
          
          % Use David Wingate's implementation:
          
          % % Blank initialization (doesn't work):
          % u = [];
          % s = [];
          % v = [];
          % First columns initialization:
          desired_rank = desired_shape_space_dimensionality + 2 + number_additional_dimensions;
          [u, s, v] = svds(distances_squared_centered(:, 1:desired_rank));
          % [u, s, v] = svd(full(distances_squared_centered(:, 1:desired_rank)), 'econ');
          for column_index = desired_rank + 1:number_data
            % whos u s v
            % [u, s, v] = addblock_svd_update(u, s, v, distances_squared_centered(:, column_index), false);
            [u, s, v] = addblock_svd_update(u, s, v, distances_squared_centered(:, column_index), true);
          end
          
        end
        
        if false
          % Debug info:
          whos u s v
          % For test "vectors_num-points00100_noise0.00_directly_embedded_method_incremental-svd" does poorly with my implementation but is close to identity for David Wingate's with either svd or svds:
          u_v_inner_products = u' * v
          beep, keyboard
        end
        
        sparse_eigenvectors = u;
        sparse_eigenvalues = diag(diag(s) .* sign(sum(u .* v, 1)).');
        
        unfiltered_sparse_eigenvalues = diag(sparse_eigenvalues).';
        embedding_info.unfiltered_sparse_eigenvalues = unfiltered_sparse_eigenvalues;
        
        good_eigenvalues = abs(unfiltered_sparse_eigenvalues) > eps(max(abs(unfiltered_sparse_eigenvalues))) * 1e2;
        
        if embedding_info.force_positive_definiteness
          unfiltered_eigenvalues_negative_proportion = abs(sum(unfiltered_sparse_eigenvalues(unfiltered_sparse_eigenvalues < 0))) / sum(abs(unfiltered_sparse_eigenvalues));
          
          good_eigenvalues = good_eigenvalues & (unfiltered_sparse_eigenvalues > 0);
          
          embedding_info.unfiltered_eigenvalues_negative_proportion = unfiltered_eigenvalues_negative_proportion;
        end

        % Respect desired_shape_space_dimensionality:
        % good_eigenvalues(desired_shape_space_dimensionality + 1:end) = false;
        good_eigenvalues(cumsum(good_eigenvalues) > desired_shape_space_dimensionality) = false;
        
        % Apply eigenvalue filtering:
        sparse_eigenvalues = diag(unfiltered_sparse_eigenvalues .* good_eigenvalues);
        
        eigenvalue_reordering = [find(good_eigenvalues), find(~good_eigenvalues)];
        sparse_eigenvalues = sparse_eigenvalues(eigenvalue_reordering, eigenvalue_reordering);
        sparse_eigenvectors = sparse_eigenvectors(:, eigenvalue_reordering);
        good_eigenvalues = good_eigenvalues(eigenvalue_reordering);
        sparse_eigenvalues = sparse_eigenvalues(good_eigenvalues, good_eigenvalues);
        sparse_eigenvectors = sparse_eigenvectors(:, good_eigenvalues);

        current_dimensionality = sum(good_eigenvalues);
        
        number_additional_dimensions = number_additional_dimensions + (desired_shape_space_dimensionality - current_dimensionality);
        
      end
        
        
      % keyboard
      % With a good SVD we should be able to detect if u and v are pointing in the same or opposite directions...
      % With a good SVD we should be able to detect if u and v are pointing in the same or opposite directions...
      positions = sparse_eigenvectors * sqrt(abs(sparse_eigenvalues));
      positions = positions(:, 1:min(end, desired_shape_space_dimensionality));
      sparse_eigenvalues = sparse_eigenvalues(1:size(positions, 2), 1:size(positions, 2));
      % mass_matrix = eye(size(positions, 2));
      mass_matrix = diag(sign(diag(sparse_eigenvalues)));
        
      
      
      
    case 'svt'
      
      
      % 2014-01-06 tebuck: Uses Luong's singular value thresholding code.
      % error('Not yet implemented!')
      
      % This assumes observations exist for every row and column and only centers observed distances:
      distances_squared_centered = distances.^2;
      row_means = sum(distances_squared_centered, 2) ./ sum(given_entries_with_diagonals, 2);
      column_means = sum(distances_squared_centered, 1) ./ sum(given_entries_with_diagonals, 1);
      overall_mean = mean(distances_squared_centered(given_entries_with_diagonals));
      [rows, columns] = find(given_entries_with_diagonals);
      for entry_index = 1:length(rows)
        row = rows(entry_index);
        column = columns(entry_index);
        % Platt 2005 (1):
        distances_squared_centered(row, column) = -1/2 * (distances_squared_centered(row, column) - row_means(row) - column_means(column) + overall_mean);
      end


      current_dimensionality = 0;
      number_additional_dimensions = 0;
      
      while current_dimensionality < desired_shape_space_dimensionality && (desired_shape_space_dimensionality + 2 + number_additional_dimensions <= number_data)
        
        % [~, ~, u, s, v] = SVT(distances_squared_centered, desired_shape_space_dimensionality + 2, 1e-4, struct('incrl', 1, 'tau_scale', 5 / number_data));
        % [~, ~, u, s, v] = SVT(distances_squared_centered, desired_shape_space_dimensionality + 2 + number_additional_dimensions, 1e-4, struct('incrl', 1, 'tau_scale', 5 / number_data));
        % % Setting tau to zero, let's see if the lack of thresholding produces acceptable results:
        % [~, ~, u, s, v] = SVT(distances_squared_centered, desired_shape_space_dimensionality + 2 + number_additional_dimensions, 1e-4, struct('incrl', 1, 'tau_scale', 0));
        % % Try a higher value of tau related to the estimate of the largest singular value of the matrix with missing entries. Assume that we are only interested in singular values with more than about 5% of the maximum:
        % [~, ~, u, s, v] = SVT(distances_squared_centered, desired_shape_space_dimensionality + 2 + number_additional_dimensions, 1e-4, struct('incrl', 1, 'tau_scale', normest(distances_squared_centered) * 5e-2 / number_data));
        % Try a higher value of tau related to the estimate of the largest singular value of the matrix with missing entries. Assume that we are only interested in singular values with more than about 5% of the maximum, and adjust the step size:
        default_step_size_scale = 1.2;
        default_step_size = default_step_size_scale * (numel(distances_squared_centered) / (numel(distances_squared_centered) - sum(given_entries_with_diagonals(:))));
        maximum_step_size = .25;
        % maximum_step_size = 1;
        % maximum_step_size = 2;
        maximum_step_size_scale = maximum_step_size / (numel(distances_squared_centered) / (numel(distances_squared_centered) - sum(given_entries_with_diagonals(:))));
        maximum_step_size = maximum_step_size_scale * (numel(distances_squared_centered) / (numel(distances_squared_centered) - sum(given_entries_with_diagonals(:))));
        % fprintf('maximum_step_size_scale = %.2e\n', maximum_step_size_scale)
        [~, ~, u, s, v] = SVT(distances_squared_centered, desired_shape_space_dimensionality + 2 + number_additional_dimensions, 1e-4, struct('incrl', 1, 'tau_scale', normest(distances_squared_centered) * 5e-2 / number_data, 'stepsize_scale', min(1.2, maximum_step_size_scale)));

        reconstructed = double(given_entries);
        [observation_rows, observation_columns] = find(given_entries);
        for entry_index = 1:length(observation_rows)
          row = observation_rows(entry_index);
          column = observation_columns(entry_index);
          reconstructed(row, column) = u(row, :) * s * v(column, :)';
        end
        reconstruction_frobenius_error_ratio = norm(reconstructed(given_entries) - distances_squared_centered(given_entries), 'fro') ./ norm(distances_squared_centered(given_entries), 'fro');
        % reconstruction_frobenius_error_ratio
        % fprintf('reconstruction_frobenius_error_ratio = %.2e\n', reconstruction_frobenius_error_ratio)
        
        sparse_eigenvectors = u;
        sparse_eigenvalues = diag(diag(s).' .* sign(sum(u .* v, 1)));
        
        unfiltered_sparse_eigenvalues = diag(sparse_eigenvalues).';
        embedding_info.unfiltered_sparse_eigenvalues = unfiltered_sparse_eigenvalues;
        
        good_eigenvalues = abs(unfiltered_sparse_eigenvalues) > eps(max(abs(unfiltered_sparse_eigenvalues))) * 1e2;
        
        if embedding_info.force_positive_definiteness
          unfiltered_eigenvalues_negative_proportion = abs(sum(unfiltered_sparse_eigenvalues(unfiltered_sparse_eigenvalues < 0))) / sum(abs(unfiltered_sparse_eigenvalues));
          
          good_eigenvalues = good_eigenvalues & (unfiltered_sparse_eigenvalues > 0);
          
          embedding_info.unfiltered_eigenvalues_negative_proportion = unfiltered_eigenvalues_negative_proportion;
        end

        % Respect desired_shape_space_dimensionality:
        % good_eigenvalues(desired_shape_space_dimensionality + 1:end) = false;
        good_eigenvalues(cumsum(good_eigenvalues) > desired_shape_space_dimensionality) = false;
        
        % Apply eigenvalue filtering:
        sparse_eigenvalues = diag(unfiltered_sparse_eigenvalues .* good_eigenvalues);
        
        eigenvalue_reordering = [find(good_eigenvalues), find(~good_eigenvalues)];
        sparse_eigenvalues = sparse_eigenvalues(eigenvalue_reordering, eigenvalue_reordering);
        sparse_eigenvectors = sparse_eigenvectors(:, eigenvalue_reordering);
        good_eigenvalues = good_eigenvalues(eigenvalue_reordering);
        sparse_eigenvalues = sparse_eigenvalues(good_eigenvalues, good_eigenvalues);
        sparse_eigenvectors = sparse_eigenvectors(:, good_eigenvalues);

        current_dimensionality = sum(good_eigenvalues);
        
        number_additional_dimensions = number_additional_dimensions + (desired_shape_space_dimensionality - current_dimensionality);
        
      end

      
      % keyboard
      % With a good SVD we should be able to detect if u and v are pointing in the same or opposite directions...
      % With a good SVD we should be able to detect if u and v are pointing in the same or opposite directions...
      positions = sparse_eigenvectors * sqrt(abs(sparse_eigenvalues));
      positions = positions(:, 1:min(end, desired_shape_space_dimensionality));
      sparse_eigenvalues = sparse_eigenvalues(1:size(positions, 2), 1:size(positions, 2));
      % mass_matrix = eye(size(positions, 2));
      mass_matrix = sign(sparse_eigenvalues);
      
      
    case 'evt'
      
      
      % 2014-01-11 tebuck: Uses Luong's singular value thresholding code modified for what I'm calling "eigenvalue thresholding" (EVT).
      % error('Not yet implemented!')
      
      % This assumes observations exist for every row and column and only centers observed distances:
      % Can we avoid this entirely by decomposing the squared distances and transforming later?
      % Can we avoid this entirely by decomposing the squared distances and transforming later?
      % Can we avoid this entirely by decomposing the squared distances and transforming later?
      % Can we avoid this entirely by decomposing the squared distances and transforming later?
      % Can we avoid this entirely by decomposing the squared distances and transforming later?
      % Can we avoid this entirely by decomposing the squared distances and transforming later?
      distances_squared_centered = distances.^2;
      row_means = sum(distances_squared_centered, 2) ./ sum(given_entries_with_diagonals, 2);
      column_means = sum(distances_squared_centered, 1) ./ sum(given_entries_with_diagonals, 1);
      overall_mean = mean(distances_squared_centered(given_entries_with_diagonals));
      [rows, columns] = find(given_entries_with_diagonals);
      for entry_index = 1:length(rows)
        row = rows(entry_index);
        column = columns(entry_index);
        % Platt 2005 (1):
        distances_squared_centered(row, column) = -1/2 * (distances_squared_centered(row, column) - row_means(row) - column_means(column) + overall_mean);
      end


      current_dimensionality = 0;
      number_additional_dimensions = 0;
      
      while current_dimensionality < desired_shape_space_dimensionality && (desired_shape_space_dimensionality + 2 + number_additional_dimensions <= number_data)
        
        % Parameter values come from experiments with method svt:
        default_step_size_scale = 1.2;
        default_step_size = default_step_size_scale * (numel(distances_squared_centered) / (numel(distances_squared_centered) - sum(given_entries_with_diagonals(:))));
        maximum_step_size = .25;
        maximum_step_size_scale = maximum_step_size / (numel(distances_squared_centered) / (numel(distances_squared_centered) - sum(given_entries_with_diagonals(:))));
        maximum_step_size = maximum_step_size_scale * (numel(distances_squared_centered) / (numel(distances_squared_centered) - sum(given_entries_with_diagonals(:))));
        % fprintf('maximum_step_size_scale = %.2e\n', maximum_step_size_scale)
        [~, ~, u, s] = SVT(distances_squared_centered, desired_shape_space_dimensionality + 2 + number_additional_dimensions, 1e-4, struct('incrl', 1, 'tau_scale', normest(distances_squared_centered) * 5e-2 / number_data, 'stepsize_scale', min(1.2, maximum_step_size_scale), 'symmetric', true));
        % % Try a threshold of zero:
        % [~, ~, u, s] = SVT(distances_squared_centered, desired_shape_space_dimensionality + 2 + number_additional_dimensions, 1e-4, struct('incrl', 1, 'tau_scale', 0, 'stepsize_scale', min(1.2, maximum_step_size_scale), 'symmetric', true));

        reconstructed = double(given_entries);
        [observation_rows, observation_columns] = find(given_entries);
        for entry_index = 1:length(observation_rows)
          row = observation_rows(entry_index);
          column = observation_columns(entry_index);
          reconstructed(row, column) = u(row, :) * s * u(column, :)';
        end
        reconstruction_frobenius_error_ratio = norm(reconstructed(given_entries) - distances_squared_centered(given_entries), 'fro') ./ norm(distances_squared_centered(given_entries), 'fro');
        % reconstruction_frobenius_error_ratio
        % fprintf('reconstruction_frobenius_error_ratio = %.2e\n', reconstruction_frobenius_error_ratio)
        
        sparse_eigenvectors = u;
        sparse_eigenvalues = s;
        
        unfiltered_sparse_eigenvalues = diag(sparse_eigenvalues).';
        embedding_info.unfiltered_sparse_eigenvalues = unfiltered_sparse_eigenvalues;
        
        good_eigenvalues = abs(unfiltered_sparse_eigenvalues) > eps(max(abs(unfiltered_sparse_eigenvalues))) * 1e2;
        
        if embedding_info.force_positive_definiteness
          unfiltered_eigenvalues_negative_proportion = abs(sum(unfiltered_sparse_eigenvalues(unfiltered_sparse_eigenvalues < 0))) / sum(abs(unfiltered_sparse_eigenvalues));
          
          good_eigenvalues = good_eigenvalues & (unfiltered_sparse_eigenvalues > 0);
          
          embedding_info.unfiltered_eigenvalues_negative_proportion = unfiltered_eigenvalues_negative_proportion;
        end

        % Respect desired_shape_space_dimensionality:
        % good_eigenvalues(desired_shape_space_dimensionality + 1:end) = false;
        good_eigenvalues(cumsum(good_eigenvalues) > desired_shape_space_dimensionality) = false;
        
        % Apply eigenvalue filtering:
        sparse_eigenvalues = diag(unfiltered_sparse_eigenvalues .* good_eigenvalues);
        
        eigenvalue_reordering = [find(good_eigenvalues), find(~good_eigenvalues)];
        sparse_eigenvalues = sparse_eigenvalues(eigenvalue_reordering, eigenvalue_reordering);
        sparse_eigenvectors = sparse_eigenvectors(:, eigenvalue_reordering);
        good_eigenvalues = good_eigenvalues(eigenvalue_reordering);
        sparse_eigenvalues = sparse_eigenvalues(good_eigenvalues, good_eigenvalues);
        sparse_eigenvectors = sparse_eigenvectors(:, good_eigenvalues);

        current_dimensionality = sum(good_eigenvalues);
        
        number_additional_dimensions = number_additional_dimensions + (desired_shape_space_dimensionality - current_dimensionality);
        
      end

      
      % keyboard
      % With a good SVD we should be able to detect if u and v are pointing in the same or opposite directions...
      % With a good SVD we should be able to detect if u and v are pointing in the same or opposite directions...
      positions = sparse_eigenvectors * sqrt(abs(sparse_eigenvalues));
      positions = positions(:, 1:min(end, desired_shape_space_dimensionality));
      sparse_eigenvalues = sparse_eigenvalues(1:size(positions, 2), 1:size(positions, 2));
      % mass_matrix = eye(size(positions, 2));
      mass_matrix = sign(sparse_eigenvalues);
      
      
    case 'eigs'
      
      
      % 2014-01-06 tebuck: Just run eigs on the sparse similarity matrix.
      % error('Not yet implemented!')
      
      % This assumes observations exist for every row and column and only centers observed distances:
      distances_squared_centered = distances.^2;
      row_means = sum(distances_squared_centered, 2) ./ sum(given_entries_with_diagonals, 2);
      column_means = sum(distances_squared_centered, 1) ./ sum(given_entries_with_diagonals, 1);
      overall_mean = mean(distances_squared_centered(given_entries_with_diagonals));
      [rows, columns] = find(given_entries_with_diagonals);
      for entry_index = 1:length(rows)
        row = rows(entry_index);
        column = columns(entry_index);
        % Platt 2005 (1):
        distances_squared_centered(row, column) = -1/2 * (distances_squared_centered(row, column) - row_means(row) - column_means(column) + overall_mean);
      end
      
      current_dimensionality = 0;
      number_additional_dimensions = 0;
      
      while current_dimensionality < desired_shape_space_dimensionality && (desired_shape_space_dimensionality + 2 + number_additional_dimensions <= number_data)
        
        % [u, d] = eigs(distances_squared_centered, desired_shape_space_dimensionality + 2);
        [sparse_eigenvectors, sparse_eigenvalues] = eigs(distances_squared_centered, desired_shape_space_dimensionality + 2 + number_additional_dimensions);

        
        unfiltered_sparse_eigenvalues = diag(sparse_eigenvalues).';
        embedding_info.unfiltered_sparse_eigenvalues = unfiltered_sparse_eigenvalues;
        
        good_eigenvalues = abs(unfiltered_sparse_eigenvalues) > eps(max(abs(unfiltered_sparse_eigenvalues))) * 1e2;
        
        if embedding_info.force_positive_definiteness
          unfiltered_eigenvalues_negative_proportion = abs(sum(unfiltered_sparse_eigenvalues(unfiltered_sparse_eigenvalues < 0))) / sum(abs(unfiltered_sparse_eigenvalues));
          
          good_eigenvalues = good_eigenvalues & (unfiltered_sparse_eigenvalues > 0);
          
          embedding_info.unfiltered_eigenvalues_negative_proportion = unfiltered_eigenvalues_negative_proportion;
        end

        % Respect desired_shape_space_dimensionality:
        % good_eigenvalues(desired_shape_space_dimensionality + 1:end) = false;
        good_eigenvalues(cumsum(good_eigenvalues) > desired_shape_space_dimensionality) = false;
        
        % Apply eigenvalue filtering:
        sparse_eigenvalues = diag(unfiltered_sparse_eigenvalues .* good_eigenvalues);
        
        eigenvalue_reordering = [find(good_eigenvalues), find(~good_eigenvalues)];
        sparse_eigenvalues = sparse_eigenvalues(eigenvalue_reordering, eigenvalue_reordering);
        sparse_eigenvectors = sparse_eigenvectors(:, eigenvalue_reordering);
        good_eigenvalues = good_eigenvalues(eigenvalue_reordering);
        sparse_eigenvalues = sparse_eigenvalues(good_eigenvalues, good_eigenvalues);
        sparse_eigenvectors = sparse_eigenvectors(:, good_eigenvalues);

        current_dimensionality = sum(good_eigenvalues);
        
        number_additional_dimensions = number_additional_dimensions + (desired_shape_space_dimensionality - current_dimensionality);
        
      end
      
      % keyboard
      positions = sparse_eigenvectors * sqrt(abs(sparse_eigenvalues));
      positions = positions(:, 1:min(end, desired_shape_space_dimensionality));
      sparse_eigenvalues = sparse_eigenvalues(1:size(positions, 2), 1:size(positions, 2));
      % mass_matrix = eye(size(positions, 2));
      mass_matrix = sign(sparse_eigenvalues);
      
      
    % case 'fpca'
      
      
      % % 2014-01-08 tebuck: New code downloaded from author for Ma et al 2011, "Fixed point and Bregman iterative methods for matrix rank minimization:"
      % error('Not yet implemented!')
      
      % % This assumes observations exist for every row and column and only centers observed distances:
      % distances_squared_centered = distances.^2;
      % row_means = sum(distances_squared_centered, 2) ./ sum(given_entries_with_diagonals, 2);
      % column_means = sum(distances_squared_centered, 1) ./ sum(given_entries_with_diagonals, 1);
      % overall_mean = mean(distances_squared_centered(given_entries_with_diagonals));
      % [rows, columns] = find(given_entries_with_diagonals);
      % for entry_index = 1:length(rows)
        % row = rows(entry_index);
        % column = columns(entry_index);
        % % Platt 2005 (1):
        % distances_squared_centered(row, column) = -1/2 * (distances_squared_centered(row, column) - row_means(row) - column_means(column) + overall_mean);
      % end
      
      % current_dimensionality = 0;
      % number_additional_dimensions = 0;
      
      % while current_dimensionality < desired_shape_space_dimensionality && (desired_shape_space_dimensionality + 2 + number_additional_dimensions <= number_data)
        
        % % % [u, d] = eigs(distances_squared_centered, desired_shape_space_dimensionality + 2);
        % % [sparse_eigenvectors, sparse_eigenvalues] = eigs(distances_squared_centered, desired_shape_space_dimensionality + 2 + number_additional_dimensions);
        
        
        % % Adapted from One_run.m from author:
        % % fr is the freedom of set of rank-r matrix, maxr is the maximum rank one
        % % can recover with p samples, which is the max rank to keep fr < 1
        % dimensionality_to_attempt = desired_shape_space_dimensionality + 2 + number_additional_dimensions;
        % p = sum(given_entries_with_diagonals(:));
        % sr = p ./ number_data^2;
        % fr = dimensionality_to_attempt*(number_data*2-dimensionality_to_attempt)/p;
        % maxr = floor(((number_data*2)-sqrt((number_data*2)^2-4*p))/2);
        % fpca_options = get_opts_FPCA(distances_squared_centered, maxr, number_data, number_data, sr, fr); 
        % fpca_result = FPCA_MatComp(number_data, number_data, find(given_entries_with_diagonals), distances_squared_centered(given_entries_with_diagonals), fpca_options);
        
        % if ~false
          % % Debug info:
          % fpca_options
          % keyboard
        % end

        
        % unfiltered_sparse_eigenvalues = diag(sparse_eigenvalues).';
        % embedding_info.unfiltered_sparse_eigenvalues = unfiltered_sparse_eigenvalues;
        
        % good_eigenvalues = abs(unfiltered_sparse_eigenvalues) > eps(max(abs(unfiltered_sparse_eigenvalues))) * 1e2;
        
        % if options.force_positive_definiteness
          % unfiltered_eigenvalues_negative_proportion = abs(sum(unfiltered_sparse_eigenvalues(unfiltered_sparse_eigenvalues < 0))) / sum(abs(unfiltered_sparse_eigenvalues));
          
          % good_eigenvalues = good_eigenvalues & (unfiltered_sparse_eigenvalues > 0);
          
          % embedding_info.unfiltered_eigenvalues_negative_proportion = unfiltered_eigenvalues_negative_proportion;
        % end

        % % Respect desired_shape_space_dimensionality:
        % % good_eigenvalues(desired_shape_space_dimensionality + 1:end) = false;
        % good_eigenvalues(cumsum(good_eigenvalues) > desired_shape_space_dimensionality) = false;
        
        % % Apply eigenvalue filtering:
        % sparse_eigenvalues = diag(unfiltered_sparse_eigenvalues .* good_eigenvalues);
        
        % eigenvalue_reordering = [find(good_eigenvalues), find(~good_eigenvalues)];
        % sparse_eigenvalues = sparse_eigenvalues(eigenvalue_reordering, eigenvalue_reordering);
        % sparse_eigenvectors = sparse_eigenvectors(:, eigenvalue_reordering);
        % good_eigenvalues = good_eigenvalues(eigenvalue_reordering);
        % sparse_eigenvalues = sparse_eigenvalues(good_eigenvalues, good_eigenvalues);
        % sparse_eigenvectors = sparse_eigenvectors(:, good_eigenvalues);

        % current_dimensionality = sum(good_eigenvalues);
        
        % number_additional_dimensions = number_additional_dimensions + (desired_shape_space_dimensionality - current_dimensionality);
        
      % end
      
      % % keyboard
      % positions = sparse_eigenvectors * sqrt(abs(sparse_eigenvalues));
      % positions = positions(:, 1:min(end, desired_shape_space_dimensionality));
      % sparse_eigenvalues = sparse_eigenvalues(1:size(positions, 2), 1:size(positions, 2));
      % % mass_matrix = eye(size(positions, 2));
      % mass_matrix = sign(sparse_eigenvalues);
      
      
    case 'lmafit'
      
      
      % 2014-01-08 tebuck: New code downloaded from author for Ma et al 2011, "Fixed point and Bregman iterative methods for matrix rank minimization:"
      % error('Not yet implemented!')
      
      % This assumes observations exist for every row and column and only centers observed distances:
      distances_squared_centered = distances.^2;
      row_means = sum(distances_squared_centered, 2) ./ sum(given_entries_with_diagonals, 2);
      column_means = sum(distances_squared_centered, 1) ./ sum(given_entries_with_diagonals, 1);
      overall_mean = mean(distances_squared_centered(given_entries_with_diagonals));
      [rows, columns] = find(given_entries_with_diagonals);
      for entry_index = 1:length(rows)
        row = rows(entry_index);
        column = columns(entry_index);
        % Platt 2005 (1):
        distances_squared_centered(row, column) = -1/2 * (distances_squared_centered(row, column) - row_means(row) - column_means(column) + overall_mean);
      end
      
      current_dimensionality = 0;
      number_additional_dimensions = 0;
      
      while current_dimensionality < desired_shape_space_dimensionality && (desired_shape_space_dimensionality + 2 + number_additional_dimensions <= number_data)
        
        % % [u, d] = eigs(distances_squared_centered, desired_shape_space_dimensionality + 2);
        % [sparse_eigenvectors, sparse_eigenvalues] = eigs(distances_squared_centered, desired_shape_space_dimensionality + 2 + number_additional_dimensions);
        
        
        % Adapted from One_run.m from author:
        % fr is the freedom of set of rank-r matrix, maxr is the maximum rank one
        % can recover with p samples, which is the max rank to keep fr < 1
        dimensionality_to_attempt = desired_shape_space_dimensionality + 2 + number_additional_dimensions;
        p = sum(given_entries_with_diagonals(:));
        sr = p ./ number_data^2;
        fr = dimensionality_to_attempt*(number_data*2-dimensionality_to_attempt)/p;
        maxr = floor(((number_data*2)-sqrt((number_data*2)^2-4*p))/2);
        % fpca_options = get_opts_FPCA(distances_squared_centered, maxr, number_data, number_data, sr, fr); 
        % fpca_result = FPCA_MatComp(number_data, number_data, find(given_entries_with_diagonals), distances_squared_centered(given_entries_with_diagonals), fpca_options);
        lmafit_options = struct;
        lmafit_options.Zfull = 0;
        lmafit_options.est_rank = 2;
        lmafit_options.rank_max = dimensionality_to_attempt;
        lmafit_options.rk_inc = 1;
        lmafit_options.print = 2;
        lmafit_options.init = 1;
        [u, s, v] = svds(distances_squared_centered, dimensionality_to_attempt);
        initial_sparse_eigenvectors = u;
        initial_sparse_eigenvalues = diag(diag(s) .* sign(sum(u .* v, 1)).');
        lmafit_options.X = u * sqrt(abs(s));
        lmafit_options.Y = (u * sqrt(abs(s)) * diag(sign(diag(s)))).';
        % [lmafit_factor1, lmafit_factor2, lmafit_result_info] = lmafit_mc_adp(number_data, number_data, dimensionality_to_attempt, find(given_entries_with_diagonals), distances_squared_centered(given_entries_with_diagonals), lmafit_options);
        [lmafit_factor1, lmafit_factor2, lmafit_result_info] = lmafit_mc_adp(number_data, number_data, 1, find(given_entries_with_diagonals), distances_squared_centered(given_entries_with_diagonals), lmafit_options);
        
        if ~false
          % Debug info:
          lmafit_options
          lmafit_result_info
          keyboard
        end

        
        unfiltered_sparse_eigenvalues = diag(sparse_eigenvalues).';
        embedding_info.unfiltered_sparse_eigenvalues = unfiltered_sparse_eigenvalues;
        
        good_eigenvalues = abs(unfiltered_sparse_eigenvalues) > eps(max(abs(unfiltered_sparse_eigenvalues))) * 1e2;
        
        if embedding_info.force_positive_definiteness
          unfiltered_eigenvalues_negative_proportion = abs(sum(unfiltered_sparse_eigenvalues(unfiltered_sparse_eigenvalues < 0))) / sum(abs(unfiltered_sparse_eigenvalues));
          
          good_eigenvalues = good_eigenvalues & (unfiltered_sparse_eigenvalues > 0);
          
          embedding_info.unfiltered_eigenvalues_negative_proportion = unfiltered_eigenvalues_negative_proportion;
        end

        % Respect desired_shape_space_dimensionality:
        % good_eigenvalues(desired_shape_space_dimensionality + 1:end) = false;
        good_eigenvalues(cumsum(good_eigenvalues) > desired_shape_space_dimensionality) = false;
        
        % Apply eigenvalue filtering:
        sparse_eigenvalues = diag(unfiltered_sparse_eigenvalues .* good_eigenvalues);
        
        eigenvalue_reordering = [find(good_eigenvalues), find(~good_eigenvalues)];
        sparse_eigenvalues = sparse_eigenvalues(eigenvalue_reordering, eigenvalue_reordering);
        sparse_eigenvectors = sparse_eigenvectors(:, eigenvalue_reordering);
        good_eigenvalues = good_eigenvalues(eigenvalue_reordering);
        sparse_eigenvalues = sparse_eigenvalues(good_eigenvalues, good_eigenvalues);
        sparse_eigenvectors = sparse_eigenvectors(:, good_eigenvalues);

        current_dimensionality = sum(good_eigenvalues);
        
        number_additional_dimensions = number_additional_dimensions + (desired_shape_space_dimensionality - current_dimensionality);
        
      end
      
      % keyboard
      positions = sparse_eigenvectors * sqrt(abs(sparse_eigenvalues));
      positions = positions(:, 1:min(end, desired_shape_space_dimensionality));
      sparse_eigenvalues = sparse_eigenvalues(1:size(positions, 2), 1:size(positions, 2));
      % mass_matrix = eye(size(positions, 2));
      mass_matrix = sign(sparse_eigenvalues);
      
      
    case 'mishra'
      
      
      [rows, columns] = find(given_entries);
      good_indices = rows > columns;
      rows = rows(good_indices);
      columns = columns(good_indices);
      % 2014-01-09 tebuck: Oops, bad indexing. :(
      % given_distances = distances(good_indices);
      given_distances = distances(sub2ind(size(given_entries), rows, columns));
      given_distances = full(given_distances);
      % Their demo code gives lowrank_dist_completion squared known distances, the code assumes it, and the paper says "A Euclidean distance matrix [...] contains the (squared) pairwise distances between n data points:"
      % Produces bad results with current options 00:35 2013-12-30...
      given_distances = given_distances.^2;
      
      % Paper suggest global solution occurs only when starting from rank 1, so do that:
      initial_positions = randn(number_data, 1);
      % initial_positions = randn(number_data, mishra_parameters.pmax);
      
      [mishra_data, info_structure] = lowrank_dist_completion(mishra_method_function, rows, columns, given_distances, initial_positions, mishra_parameters);

      positions = mishra_data(:, 1:min(end, desired_shape_space_dimensionality));
      mass_matrix = eye(size(positions, 2));
      
      
    case 'mishra-svt-init'
      
      
      [rows, columns] = find(given_entries);
      good_indices = rows > columns;
      rows = rows(good_indices);
      columns = columns(good_indices);
      % 2014-01-09 tebuck: Oops, bad indexing. :(
      % given_distances = distances(good_indices);
      given_distances = distances(sub2ind(size(given_entries), rows, columns));
      given_distances = full(given_distances);
      % Their demo code gives lowrank_dist_completion squared known distances, the code assumes it, and the paper says "A Euclidean distance matrix [...] contains the (squared) pairwise distances between n data points:"
      % Produces bad results with current options 00:35 2013-12-30...
      given_distances = given_distances.^2;
      
      % % Paper suggest global solution occurs only when starting from rank 1, but initialize with SVT to see if that gives good results:
      initial_positions = embed_partial_distance_matrix(distances, setfield(embedding_info, 'method', 'svt'));
      
      [mishra_data, info_structure] = lowrank_dist_completion(mishra_method_function, rows, columns, given_distances, initial_positions, mishra_parameters);

      positions = mishra_data(:, 1:min(end, desired_shape_space_dimensionality));
      mass_matrix = eye(size(positions, 2));
      
      
    case 'mishra-eigs-init'
      
      
      [rows, columns] = find(given_entries);
      good_indices = rows > columns;
      rows = rows(good_indices);
      columns = columns(good_indices);
      % 2014-01-09 tebuck: Oops, bad indexing. :(
      % given_distances = distances(good_indices);
      given_distances = distances(sub2ind(size(given_entries), rows, columns));
      given_distances = full(given_distances);
      % Their demo code gives lowrank_dist_completion squared known distances, the code assumes it, and the paper says "A Euclidean distance matrix [...] contains the (squared) pairwise distances between n data points:"
      % Produces bad results with current options 00:35 2013-12-30...
      given_distances = given_distances.^2;
      
      % % Paper suggest global solution occurs only when starting from rank 1, but initialize with SVT to see if that gives good results:
      initial_positions = embed_partial_distance_matrix(distances, setfield(embedding_info, 'method', 'eigs'));
      
      [mishra_data, info_structure] = lowrank_dist_completion(mishra_method_function, rows, columns, given_distances, initial_positions, mishra_parameters);

      positions = mishra_data(:, 1:min(end, desired_shape_space_dimensionality));
      mass_matrix = eye(size(positions, 2));
      
      
    case 'mishra-incremental-svd-init'
      
      
      [rows, columns] = find(given_entries);
      good_indices = rows > columns;
      rows = rows(good_indices);
      columns = columns(good_indices);
      % 2014-01-09 tebuck: Oops, bad indexing. :(
      % given_distances = distances(good_indices);
      given_distances = distances(sub2ind(size(given_entries), rows, columns));
      given_distances = full(given_distances);
      % Their demo code gives lowrank_dist_completion squared known distances, the code assumes it, and the paper says "A Euclidean distance matrix [...] contains the (squared) pairwise distances between n data points:"
      % Produces bad results with current options 00:35 2013-12-30...
      given_distances = given_distances.^2;
      
      % % Paper suggest global solution occurs only when starting from rank 1, but initialize with SVT to see if that gives good results:
      initial_positions = embed_partial_distance_matrix(distances, setfield(embedding_info, 'method', 'incremental-svd'));
      
      [mishra_data, info_structure] = lowrank_dist_completion(mishra_method_function, rows, columns, given_distances, initial_positions, mishra_parameters);

      positions = mishra_data(:, 1:min(end, desired_shape_space_dimensionality));
      mass_matrix = eye(size(positions, 2));
      
      
    case 'nystrom-euclidean'
      
      
      % error('Not yet implemented!')
      % warning('Should respect options.desired_shape_space_dimensionality!')
      
      % Use only landmarks in the Nystrom approximation (see Platt 2005, "FastMap, MetricMap, and Landmark MDS are all Nystrom Algorithms"):
      
      % if ~isempty(options.distances_to_compute)
        % error('method "%s" does not support option distances_to_compute being nonempty', options.method)
      % end
      
      % Obtain A and B (variable names from Platt 2005):
      landmark_distances = distances(landmark_indices, landmark_indices);
      remainder_distances = distances(~landmark_indices, landmark_indices);
      
      % Center A and B as in Landmark MDS:
      % A, Platt 2005, equation (13):
      landmark_centering_matrix = eye(number_landmarks) - 1 / number_landmarks;
      % This becomes landmark_centering_matrix * [E_ij^2 - 1/m * sum_q E_iq^2] = [(E_ij^2 - 1/m * sum_q E_iq^2) - sum_p (E_pj^2 - 1/m * sum_q E_pq^2)]?
      landmark_kernel = -1 / 2 * landmark_centering_matrix * (landmark_distances.^2) * landmark_centering_matrix;
      % B, Platt 2005, equation (15):
      remainder_kernel = -1 / 2 * (remainder_distances.^2 - repmat(mean(landmark_distances.^2, 1), [number_data - number_landmarks, 1]));
      
      approximate_kernel = distances;
      approximate_kernel(landmark_indices, landmark_indices) = landmark_kernel;
      approximate_kernel(~landmark_indices, landmark_indices) = remainder_kernel;
      approximate_kernel(landmark_indices, ~landmark_indices) = remainder_kernel';
      
      [extended_kernel_eigenvectors, extended_kernel_eigenvalues] = nystrom_extension(approximate_kernel, struct('method', 'classical', 'force_positive_definiteness', embedding_info.force_positive_definiteness));
      [extended_kernel_eigenvectors, extended_kernel_eigenvalues] = filter_spectral_decomposition(extended_kernel_eigenvectors, extended_kernel_eigenvalues, false);
      positions = extended_kernel_eigenvectors * sqrt(abs(extended_kernel_eigenvalues));
      mass_matrix = sign(extended_kernel_eigenvalues);
      
      if false
        % Debug info:
        mass_matrix_diagonal = diag(mass_matrix).'
        % set(gcf, 'visible', 'on'), imagesc([imag(positions)]), colorbar
        % pause
        set(gcf, 'visible', 'on')
        scatter(positions(:, 1), positions(:, 2))
        axis image
        title(sprintf('In %s', mfilename), 'Interpreter', 'none')
        beep, dbstack
        % keyboard
        pause
      end

      
    case 'nystrom-dissimilarity-clip'
      error('Not yet implemented!')
      warning('Should respect options.desired_shape_space_dimensionality!')
      
      
    case 'nystrom-dissimilarity-orthogonalize'
    
    
      % error('Not yet implemented!')
      warning('Should respect options.desired_shape_space_dimensionality!')
    
    
      % Based on method 'nystrom-euclidean' but using orthogonalization from Belongie et al. 2002, "Spectral Partitioning with Indefinite Kernels Using the Nystrom Extension:"
      
      % Use only landmarks in the Nystrom approximation (see Platt 2005, "FastMap, MetricMap, and Landmark MDS are all Nystrom Algorithms"):
      
      % Obtain A and B (variable names from Platt 2005):
      landmark_distances = distances(landmark_indices, landmark_indices);
      remainder_distances = distances(~landmark_indices, landmark_indices);
      
      % Center A and B as in Landmark MDS:
      % A, Platt 2005, equation (13):
      landmark_centering_matrix = eye(number_landmarks) - 1 / number_landmarks;
      % This becomes landmark_centering_matrix * [E_ij^2 - 1/m * sum_q E_iq^2] = [(E_ij^2 - 1/m * sum_q E_iq^2) - sum_p (E_pj^2 - 1/m * sum_q E_pq^2)]?
      landmark_kernel = -1 / 2 * landmark_centering_matrix * (landmark_distances.^2) * landmark_centering_matrix;
      % B, Platt 2005, equation (13):
      % remainder_kernel = -1 / 2 * (remainder_distances.^2 - repmat(mean(landmark_distances.^2, 2).', [number_data - number_landmarks, 1]));
      remainder_kernel = -1 / 2 * (remainder_distances.^2 - repmat(mean(landmark_distances.^2, 1), [number_data - number_landmarks, 1]));
      
      
      % Decompose A, Platt 2005, equation (7):
      % if options.use_schur_decomposition
        % [landmark_kernel_eigenvectors, landmark_kernel_eigenvalues] = schur(landmark_kernel);
      % else
        % [landmark_kernel_eigenvectors, landmark_kernel_eigenvalues] = eig(landmark_kernel);
        [landmark_kernel_eigenvectors, landmark_kernel_eigenvalues] = eig(full(landmark_kernel));
        % [landmark_kernel_eigenvectors, landmark_kernel_eigenvalues] = eigs(landmark_kernel, size(landmark_kernel, 1));
      % end
      % keyboard

      
      % Correct eigenvalues:
      if embedding_info.force_positive_definiteness
        % Clip:
        minimum_landmark_kernel_eigenvalue = min(diag(landmark_kernel_eigenvalues))
        large_negative_eigenvalue_threshold = -mean(abs(diag(landmark_kernel_eigenvalues))) * eps * 1e2
        if any(diag(landmark_kernel_eigenvalues) < large_negative_eigenvalue_threshold)
          warning('Clipping large negative eigenvalues!')
        end
        landmark_kernel_eigenvalues(landmark_kernel_eigenvalues < 0) = 0;
      end
      
      
      % Final formula for W hat, Belongie et al. 2002, section 4:
      % projected_columns_kernel = [landmark_kernel_eigenvectors; remainder_kernel * landmark_kernel_eigenvectors * pinv(landmark_kernel_eigenvalues)];
      projected_columns_kernel = zeros(number_data, size(landmark_kernel_eigenvectors, 2));
      projected_columns_kernel(landmark_indices, :) = landmark_kernel_eigenvectors;
      projected_columns_kernel(~landmark_indices, :) = remainder_kernel * landmark_kernel_eigenvectors * pinv(landmark_kernel_eigenvalues);
      positions = projected_columns_kernel * sqrt(landmark_kernel_eigenvalues);
      % [positions_eigenvectors, positions_eigenvalues] = schur(positions.' * positions);
      % complete_kernel_eigenvectors = positions * positions_eigenvectors * sqrt(pinv(positions_eigenvalues));
      % complete_kernel = complete_kernel_eigenvectors * positions_eigenvalues * complete_kernel_eigenvectors.';
      % % Make symmetric:
      % complete_kernel = 1 / 2 * (complete_kernel + complete_kernel.');

      
      % % This should work in general according to Schleif and Gisbrecht 2013 (right before (5)):
      % complete_distances = sqrt(repmat(diag(complete_kernel), [1, number_data]) + repmat(diag(complete_kernel), [1, number_data]).' - 2 * complete_kernel);
      % % warning('>>>> HACK, setting imaginary distances to zero'), complete_distances = real(complete_distances);
      % distances = complete_distances;


      positions = positions(:, 1:min(end, desired_shape_space_dimensionality));
      mass_matrix = eye(size(positions, 2));
      
      
      % keyboard
      % error
      
      
    case 'sdm-nystrom'
    
    
      % error('Not yet implemented!')
      % error('Not yet finished!')
    
      % Based on method 'nystrom-euclidean' but uses the Nystrom extension on the squared distance matrix and converts that to a kernel by double centering. Thus this is my novel approach inspired by Belongie et al. 2002 (but not using its general reorthogonalization method yet):
      
      % Variable names from "Taraz Shape space estimation 006.docx":
      
      % Use only landmarks:
      
      % E:
      landmark_squared_distances = full(distances(landmark_indices, landmark_indices)).^2;
      % F (note same orientation as usual Nystrom extension and as in that document, different from other methods in this function):
      remainder_squared_distances = full(distances(landmark_indices, ~landmark_indices)).^2;
      % R, Upsilon:
      [landmark_eigenvectors, landmark_eigenvalues] = spectral_decomposition(landmark_squared_distances);
      
      % Z without centering or scaling:
      % kernel_factor_nonorthogonal_uncentered_unscaled = [landmark_eigenvectors; remainder_squared_distances' * landmark_eigenvectors * pinv(landmark_eigenvalues)];
      kernel_factor_nonorthogonal_uncentered_unscaled = [landmark_eigenvectors; remainder_squared_distances' * landmark_eigenvectors * sparse(pinv(full(landmark_eigenvalues)))];
      
      % Z without centering:
      kernel_factor_nonorthogonal_unscaled = (kernel_factor_nonorthogonal_uncentered_unscaled - repmat(mean(kernel_factor_nonorthogonal_uncentered_unscaled, 1), [number_data, 1]));
      % Z:
      kernel_factor_nonorthogonal = kernel_factor_nonorthogonal_unscaled * sqrt(-1/2 * landmark_eigenvalues);
      
      kernel_eigenvectors = kernel_factor_nonorthogonal_unscaled;
      kernel_eigenvalues = abs(-1/2 * landmark_eigenvalues);
      % Only filter the decomposition after it is complete:
      [kernel_eigenvectors, kernel_eigenvalues] = filter_spectral_decomposition(kernel_eigenvectors, kernel_eigenvalues, false);
      positions = kernel_eigenvectors * sqrt(kernel_eigenvalues);
      positions = positions(:, 1:min(end, desired_shape_space_dimensionality));
      mass_matrix = sign(kernel_eigenvalues);
      mass_matrix = mass_matrix(1:size(positions, 2), 1:size(positions, 2));
      
      if false
        % Debug info:
        mass_matrix_diagonal = diag(mass_matrix).'
        % set(gcf, 'visible', 'on'), imagesc([imag(positions)]), colorbar
        % pause
        set(gcf, 'visible', 'on')
        scatter(positions(:, 1), positions(:, 2))
        axis image
        title(sprintf('In %s', mfilename), 'Interpreter', 'none')
        beep, dbstack
        % keyboard
        pause
      end
      
      
    case 'sdm-nystrom-kernel-orthogonalize'
    
    
      % error('Not yet implemented!')
      % error('Not yet finished!')
    
      % Based on method 'nystrom-dissimilarity-orthogonalize' but uses the Nystrom extension on the squared distance matrix, converts that to a kernel by double centering, and orthogonalizes the kernel as in Belongie et al. 2002, "Spectral Partitioning with Indefinite Kernels Using the Nystrom Extension." Thus this is my novel approach inspired by Belongie et al. 2002 and Liu 2009, "Spectral mesh segmentation" (5.4):
      
      % Variable names from "Taraz Shape space estimation 007.docx":
      
      % Use only landmarks:
      
      if false
        % E:
        landmark_squared_distances = full(distances(landmark_indices, landmark_indices)).^2;
        % F (note same orientation as usual Nystrom extension and as in that document, different from other methods in this function):
        remainder_squared_distances = full(distances(landmark_indices, ~landmark_indices)).^2;
        
        % R, Upsilon:
        [landmark_eigenvectors, landmark_eigenvalues] = spectral_decomposition(landmark_squared_distances);
        
        % R bar without centering:
        % kernel_factor_nonorthogonal_unscaled_uncentered = [landmark_eigenvectors; remainder_squared_distances' * landmark_eigenvectors * sparse(pinv(full(landmark_eigenvalues)))];
        % error('Reorder data!!!!!!!')
        kernel_factor_nonorthogonal_unscaled_uncentered = zeros(number_data, number_landmarks);
        kernel_factor_nonorthogonal_unscaled_uncentered(landmark_indices, :) = landmark_eigenvectors;
        kernel_factor_nonorthogonal_unscaled_uncentered(~landmark_indices, :) = remainder_squared_distances' * landmark_eigenvectors * sparse(pinv(full(landmark_eigenvalues)));
      else
        squared_distances = distances.^2;
        % whos squared_distances, keyboard
        % R bar without centering, Upsilon:
        [kernel_factor_nonorthogonal_unscaled_uncentered, landmark_eigenvalues] = nystrom_extension(squared_distances, struct('method', 'classical', 'force_positive_definiteness', embedding_info.force_positive_definiteness));
      end
      % R bar:
      kernel_factor_nonorthogonal_unscaled = kernel_factor_nonorthogonal_unscaled_uncentered - repmat(mean(kernel_factor_nonorthogonal_unscaled_uncentered, 1), [number_data, 1]);
      % Upsilon bar:
      kernel_landmark_eigenvalues = -1/2 * landmark_eigenvalues;
      
      % Z:
      kernel_factor_nonorthogonal = kernel_factor_nonorthogonal_unscaled * sqrt(kernel_landmark_eigenvalues);
      
      % H, Sigma:
      [inner_product_eigenvectors, inner_product_eigenvalues] = spectral_decomposition(kernel_factor_nonorthogonal' * kernel_factor_nonorthogonal);
      
      % [inner_product_eigenvectors, inner_product_eigenvalues] = filter_spectral_decomposition(inner_product_eigenvectors, inner_product_eigenvalues, false);
      
      % V, Sigma:
      kernel_eigenvectors = kernel_factor_nonorthogonal * inner_product_eigenvectors * pinv(general_square_root(inner_product_eigenvalues));
      kernel_eigenvalues = inner_product_eigenvalues;
      % Only filter the decomposition after it is complete:
      [kernel_eigenvectors, kernel_eigenvalues] = filter_spectral_decomposition(kernel_eigenvectors, kernel_eigenvalues, false);
      
      % V * Sigma^(1/2):
      % % positions = kernel_factor_nonorthogonal * inner_product_eigenvectors * pinv(general_square_root(inner_product_eigenvalues)) * general_square_root(inner_product_eigenvalues);
      % positions = kernel_factor_nonorthogonal * inner_product_eigenvectors;
      % mass_matrix = diag(sign(diag(inner_product_eigenvalues)));
      positions = kernel_eigenvectors * general_square_root(kernel_eigenvalues);
      mass_matrix = diag(sign(diag(kernel_eigenvalues)));
      
      positions_real = real(positions);
      if any((positions_real - positions) > eps(max(abs(positions(:)))) * 1e5)
        warning('Removing large imaginary components of returned positions!')
      end
      positions = positions_real;
      
      if false
        % Debug info:
        mass_matrix_diagonal = diag(mass_matrix).'
        % set(gcf, 'visible', 'on'), imagesc([imag(positions)]), colorbar
        % pause
        set(gcf, 'visible', 'on')
        scatter(positions(:, 1), positions(:, 2))
        axis image
        title(sprintf('In %s', mfilename), 'Interpreter', 'none')
        beep, dbstack
        keyboard
        % pause
      end
      
      
    case 'approximate-svd'
    
    
      % error('Not yet implemented!')
      error('Not yet finished!')
      
      
        
      
      % Based on Nemtsov et al. 2013, "Matrix Compression using the Nystrom Method," section 2.2.2, and Platt 2005, section 4.2:
      
      % Based on code for nystrom-euclidean:
      
      % Obtain A_M and B_M (variable names from Nemtsov, note same orientation as usual Nystrom extension and as in that document, different from other methods in this function):
      landmark_squared_distances = full(distances(landmark_indices, landmark_indices).^2);
      remainder_squared_distances = full(distances(landmark_indices, ~landmark_indices).^2);
      data_order_restoration = zeros(number_data, 1);
      data_order_restoration(landmark_indices) = (1:number_landmarks)';
      data_order_restoration(~landmark_indices) = (number_landmarks + 1:number_data)';

      
      % Center A_M and B_M, equation (2.1) to produce kernel matrices as in Landmark MDS:
      
      % A_M (using Platt, equation (13)):
      landmark_squared_distances_column_means = mean(landmark_squared_distances, 1);
      landmark_squared_distances_mean = mean(landmark_squared_distances(:));
      landmark_kernel = -1 / 2 * (landmark_squared_distances - repmat(landmark_squared_distances_column_means, [number_landmarks, 1]) - repmat(landmark_squared_distances_column_means.', [1, number_landmarks]) + landmark_squared_distances_mean);
      % B_M (using Platt, equation (15)):
      remainder_kernel = -1 / 2 * (remainder_squared_distances - repmat(landmark_squared_distances_column_means.', [1, number_data - number_landmarks]));
      
      % Get extended left eigenvectors Uhat from U and Utilde, Nemtsov et al. 2013, equation (2.4), (2.6), (2.9):
      % Scale the eigenvectors as G_U, G_V:
      kernel_unscaled_left_eigenvectors = [landmark_kernel; remainder_kernel'];
      kernel_unscaled_left_eigenvectors = kernel_unscaled_left_eigenvectors(data_order_restoration, :);
      landmark_kernel_square_root = general_square_root(landmark_kernel);
      kernel_scaled_left_eigenvectors = kernel_unscaled_left_eigenvectors * pinv(landmark_kernel_square_root);
      % Compute the product's spectral decomposition F, Sigma:
      [kernel_scaled_eigenvector_product_eigenvectors, kernel_scaled_eigenvector_product_eigenvalues] = spectral_decomposition(kernel_scaled_left_eigenvectors' * kernel_scaled_left_eigenvectors, ~embedding_info.force_positive_definiteness);
      % Compute the final spectral decomposition of the Nystrom approximation:
      kernel_corrected_eigenvalues = kernel_scaled_eigenvector_product_eigenvalues;
      kernel_corrected_left_eigenvectors = kernel_scaled_left_eigenvectors * kernel_scaled_eigenvector_product_eigenvectors * pinv(sqrt(kernel_scaled_eigenvector_product_eigenvalues));
      % Imaginary values were close to zero in a C2C12 test:
      kernel_corrected_left_eigenvectors = real(kernel_corrected_left_eigenvectors);
      
      
      if ~false
        fprintf('landmark_kernel reconstruction relative Frobenius error: %e\n', norm(landmark_kernel - kernel_corrected_left_eigenvectors(landmark_indices, :) * kernel_corrected_eigenvalues * kernel_corrected_left_eigenvectors(landmark_indices, :)', 'fro') ./ norm(landmark_kernel, 'fro'))
        fprintf('kernel_corrected_left_eigenvectors orthogonality relative Frobenius error: %e\n', norm(eye(size(kernel_corrected_left_eigenvectors, 2)) - kernel_corrected_left_eigenvectors' * kernel_corrected_left_eigenvectors, 'fro') ./ norm(eye(size(kernel_corrected_left_eigenvectors, 2)), 'fro'))
        fprintf('remainder_squared_distances reconstruction relative Frobenius error: %e\n', norm(remainder_kernel - kernel_corrected_left_eigenvectors(landmark_indices, :) * kernel_corrected_eigenvalues * kernel_corrected_left_eigenvectors(~landmark_indices, :)', 'fro') ./ norm(remainder_squared_distances, 'fro'))
        kernel_corrected_left_eigenvectors
        kernel_corrected_eigenvalues
        beep, keyboard
        % pause
      end
      
      
      % Remove unnecessary dimensions:
      
      unfiltered_kernel_eigenvalues = diag(kernel_corrected_eigenvalues).';
      embedding_info.unfiltered_landmark_kernel_eigenvalues = unfiltered_kernel_eigenvalues;
      
      good_eigenvalues = abs(unfiltered_kernel_eigenvalues) > eps(max(abs(unfiltered_kernel_eigenvalues))) * 1e2;
      
      if embedding_info.force_positive_definiteness
        unfiltered_eigenvalues_negative_proportion = abs(sum(unfiltered_kernel_eigenvalues(unfiltered_kernel_eigenvalues < 0))) / sum(abs(unfiltered_kernel_eigenvalues));
        % fprintf('Removing negative eigenvalues (%.3f of "variance")\n', unfiltered_eigenvalues_negative_proportion)
        good_eigenvalues = good_eigenvalues & (unfiltered_kernel_eigenvalues > 0);
        embedding_info.unfiltered_eigenvalues_negative_proportion = unfiltered_eigenvalues_negative_proportion;
      end
      
      % Respect desired_shape_space_dimensionality:
      % good_eigenvalues(desired_shape_space_dimensionality + 1:end) = false;
      good_eigenvalues(cumsum(good_eigenvalues) > desired_shape_space_dimensionality) = false;
      
      % Apply eigenvalue filtering:
      kernel_corrected_left_eigenvectors = diag(unfiltered_kernel_eigenvalues .* good_eigenvalues);
      
      eigenvalue_reordering = [find(good_eigenvalues), find(~good_eigenvalues)];
      kernel_corrected_left_eigenvectors = kernel_corrected_left_eigenvectors(eigenvalue_reordering, eigenvalue_reordering);
      kernel_corrected_left_eigenvectors = kernel_corrected_left_eigenvectors(:, eigenvalue_reordering);
      good_eigenvalues = good_eigenvalues(eigenvalue_reordering);
      kernel_corrected_left_eigenvectors = kernel_corrected_left_eigenvectors(good_eigenvalues, good_eigenvalues);
      kernel_corrected_left_eigenvectors = kernel_corrected_left_eigenvectors(:, good_eigenvalues);

      % S hat star, Schleif and Gisbrecht 2013, equation (4) (could use equation (5) at some point to work on distances directly?):
      % columns_kernel = [landmark_kernel; remainder_kernel];
      % positions = columns_kernel * kernel_corrected_left_eigenvectors * sqrt(pinv(kernel_corrected_left_eigenvectors));
      positions = kernel_corrected_left_eigenvectors * sqrt(pinv(kernel_corrected_left_eigenvectors));
      
      mass_matrix = eye(size(positions, 2));
      
      if false
        % Debug info:
        if size(positions, 2) < 2
          beep
          keyboard
        end
      end
      
      
      if false
        % Debug info:
        visible_state = get(gcf, 'Visible');
        set(gcf, 'Visible', 'on');
        scatter(positions(:, 1), positions(:, 2))
        axis image
        title(sprintf('In %s', mfilename), 'Interpreter', 'none')
        beep
        % keyboard
        pause
        set(gcf, 'Visible', visible_state);
      end
      
      
      
    case 'pekalska-pseudo-euclidean'
    
    
      % warning('Should respect options.desired_shape_space_dimensionality!')

      % Based on Pseudo-Euclidean (Krein space?) method from Pekalska and Duin 2008, "Beyond Traditional Kernels: Classification in Two Dissimilarity-Based Representation Spaces," Section III.B:
      
      
      % Compute Krein space embedding for landmarks:
      
      
      
      % Compute indefinite Gram matrix G for representatives/landmarks:
      
      landmark_distances = distances(landmark_indices, landmark_indices);
      % J:
      landmark_centering_matrix = eye(number_landmarks) - 1 / number_landmarks;
      % G:
      landmark_kernel = -(1 / 2) * landmark_centering_matrix * (landmark_distances.^2) * landmark_centering_matrix;
      % Make symmetric (due to numerical error):
      landmark_kernel = (landmark_kernel + landmark_kernel.') ./ 2;
      
      % Decompose G:
      
      % if options.use_schur_decomposition
        % [landmark_kernel_eigenvectors, landmark_kernel_eigenvalues] = schur(landmark_kernel);
      % else
        % [landmark_kernel_eigenvectors, landmark_kernel_eigenvalues] = eig(landmark_kernel);
        [landmark_kernel_eigenvectors, landmark_kernel_eigenvalues] = eig(full(landmark_kernel));
        % [landmark_kernel_eigenvectors, landmark_kernel_eigenvalues] = eigs(landmark_kernel, size(landmark_kernel, 1));
      % end
      % keyboard
      % Sort by decreasing absolute eigenvalue:
      [~, landmark_kernel_eigenvalue_order] = sort(abs(diag(landmark_kernel_eigenvalues)), 1, 'descend');
      landmark_kernel_eigenvalues = landmark_kernel_eigenvalues(landmark_kernel_eigenvalue_order, landmark_kernel_eigenvalue_order);
      % landmark_kernel_eigenvectors = landmark_kernel_eigenvectors(landmark_kernel_eigenvalue_order, landmark_kernel_eigenvalue_order);
      landmark_kernel_eigenvectors = landmark_kernel_eigenvectors(:, landmark_kernel_eigenvalue_order);
      
      % if options.use_schur_decomposition
        % if false
          % landmark_kernel_eigenvalues
          % warning('schur sometimes returns non-diagonal eigenvalues!')
        % end
        % landmark_kernel_eigenvalues = diag(diag(landmark_kernel_eigenvalues));
        % if false
          % landmark_kernel_eigenvalues
        % end
      % end
      
      % The last eigenvalue is always zero:
      landmark_kernel_eigenvalues(end, end) = 0;
      
      
      % Remove unnecessary dimensions:
      
      
      unfiltered_landmark_kernel_eigenvalues = diag(landmark_kernel_eigenvalues).';
      embedding_info.unfiltered_landmark_kernel_eigenvalues = unfiltered_landmark_kernel_eigenvalues;
      
      if false
        warning('>>>> HACK, not removing any eigenvalues!')
        good_eigenvalues = true(size(unfiltered_landmark_kernel_eigenvalues));
      else
        % % Zero eigenvalues:
        % good_eigenvalues = diag(landmark_kernel_eigenvalues) ~= 0;
        % Small absolute value eigenvalues:
        % good_eigenvalues = abs(diag(landmark_kernel_eigenvalues)) > eps(max(abs(diag(landmark_kernel_eigenvalues)))) * 1e2;
        good_eigenvalues = abs(unfiltered_landmark_kernel_eigenvalues) > eps(max(abs(unfiltered_landmark_kernel_eigenvalues))) * 1e2;
        
        if embedding_info.force_positive_definiteness
          % landmark_kernel_eigenvalues
          % unfiltered_landmark_kernel_eigenvalues
          unfiltered_eigenvalues_negative_proportion = abs(sum(unfiltered_landmark_kernel_eigenvalues(unfiltered_landmark_kernel_eigenvalues < 0))) / sum(abs(unfiltered_landmark_kernel_eigenvalues));
          
          % fprintf('Removing negative eigenvalues (%.3f of "variance")\n', unfiltered_eigenvalues_negative_proportion)
          good_eigenvalues = good_eigenvalues & (unfiltered_landmark_kernel_eigenvalues > 0);
          
          embedding_info.unfiltered_eigenvalues_negative_proportion = unfiltered_eigenvalues_negative_proportion;
        end
      end
      
      % Respect desired_shape_space_dimensionality:
      good_eigenvalues(desired_shape_space_dimensionality + 1:end) = false;
      
      % Apply eigenvalue filtering:
      landmark_kernel_eigenvalues = diag(unfiltered_landmark_kernel_eigenvalues .* good_eigenvalues);
      
      eigenvalue_reordering = [find(good_eigenvalues), find(~good_eigenvalues)];
      landmark_kernel_eigenvalues = landmark_kernel_eigenvalues(eigenvalue_reordering, eigenvalue_reordering);
      landmark_kernel_eigenvectors = landmark_kernel_eigenvectors(:, eigenvalue_reordering);
      good_eigenvalues = good_eigenvalues(eigenvalue_reordering);
      landmark_kernel_eigenvalues = landmark_kernel_eigenvalues(good_eigenvalues, good_eigenvalues);
      landmark_kernel_eigenvectors = landmark_kernel_eigenvectors(:, good_eigenvalues);
      
      
      embedding_info.landmark_kernel_eigenvalues = landmark_kernel_eigenvalues;
      embedding_info.landmark_kernel_eigenvectors = landmark_kernel_eigenvectors;
      
      % Not used for this algorithm!
      
      % Embed landmarks:
      
      landmark_positions = landmark_kernel_eigenvectors * abs(landmark_kernel_eigenvalues).^(1/2);
      % I_pq:
      mass_matrix = sign(landmark_kernel_eigenvalues);
      
      
      % Compute projections for rest of data:
      
      
      % Compute indefinite cross-Gram matrix G_new for the rest of the data:
      remainder_distances = distances(~landmark_indices, landmark_indices);
      % G_new:
      remainder_kernel = -(1 / 2) * (remainder_distances.^2 - (1 / number_landmarks) * ones(number_data - number_landmarks, 1) * ones(1, number_landmarks) * landmark_distances.^2) * landmark_centering_matrix;
      
      % Embed remaining data:
      remainder_positions = remainder_kernel * landmark_positions * pinv(abs(landmark_kernel_eigenvalues)) * mass_matrix;
      
      
      % Final positions:
      positions = zeros(number_data, size(landmark_positions, 2));
      positions(landmark_indices, :) = landmark_positions;
      positions(~landmark_indices, :) = remainder_positions;
      
      % positions = positions(:, good_eigenvalues);
      % mass_matrix = mass_matrix(good_eigenvalues, good_eigenvalues);
      
      
      % Could we visualize these positions using the reciprocals of the negative eigenvalue/imaginary coordinates? Larger values lead to smaller inner products, so...
      
      % error('Not yet implemented!')
      
      
      if false
        % Debug info
        whos
        imagesc(imag(complete_distances)), axis image, title('imag(complete_distances)', 'Interpreter', 'none'), colormap('copper'), colorbar
        pause
        % whos, beep, keyboard
      end
      
      % if ~false
        dimensionality_before_forcing = size(positions, 2);
        if embedding_info.force_positive_definiteness
          euclidean_dimensions = sign(diag(mass_matrix)) > 0;
          positions = positions(:, euclidean_dimensions);
          mass_matrix = eye(size(positions, 2));
        end
        dimensionality_after_forcing = size(positions, 2);
      % end

      
      % keyboard
      % error
    
      
    case 'pekalska-dissimilarity'
    
    
      % warning('Should respect options.desired_shape_space_dimensionality!')

      % Based on dissimilarity space method from Pekalska and Duin 2008, "Beyond Traditional Kernels: Classification in Two Dissimilarity-Based Representation Spaces," Section III.B:
      
      
      % Compute Krein space embedding for landmarks:
      
      
      
      % Compute indefinite Gram matrix G for representatives/landmarks:
      % positions = distances(:, landmark_indices);
      positions = full(distances(:, landmark_indices));
      % use_pca = false;
      use_pca = true;
      embedding_info.use_pca = use_pca;
      if use_pca
        % PCA this:
        [principal_components, positions, principal_component_eigenvalues] = princomp(zscore(positions));
        embedding_info.use_pca = use_pca;
        embedding_info.principal_components = principal_components;
        embedding_info.principal_component_eigenvalues = principal_component_eigenvalues;
      end
      
      positions = positions(:, 1:min(end, desired_shape_space_dimensionality));
      mass_matrix = eye(size(positions, 2));
      
      % whos positions mass_matrix

      
      if false
        % Debug info
        whos
        imagesc(imag(complete_distances)), axis image, title('imag(complete_distances)', 'Interpreter', 'none'), colormap('copper'), colorbar
        pause
        % whos, beep, keyboard
      end

      
      % keyboard
      % error
    
      
    
      
    case 'iterated-soft-thresholding'
    
    
      warning('Should respect options.desired_shape_space_dimensionality!')

      % Based on ???:
      

      % Code adapted from Demo.m from <http://www.mathworks.com/matlabcentral/fileexchange/26395-matrix-completion-via-thresholding>:
      
      % T = randperm(prod(sizeX));
      % IDX = T(1:round(0.4*prod(sizeX))); % 40% sampling

      % Creating operator for selecting entries at the chosen random locations 
      % Requires Sparco Toolbox 
      % http://www.cs.ubc.ca/labs/scl/sparco/
      % M = opRestriction(prod(sizeX), IDX);
      % M = opRestriction(prod(size(given_entries)), find(given_entries(:)));
      M = opMask(given_entries | speye(size(given_entries)));

      % Sampled data
      y = M(distances(:), 1);
      whos distances M y
      completed_distances = IST_MC(y, M, size(given_entries));
      completed_distances = (completed_distances + completed_distances.') ./ 2;
      completed_distances(speye(size(given_entries)) == 1) = 0;
      completed_distances(completed_distances < 0) = 0;
      
      % End code adapted from Demo.m.

      
      positions = cmdscale(completed_distances);
      mass_matrix = eye(size(positions, 2));

      
      % keyboard
      % error
    
      
    case 'similarity-principal-components'
      
      
      % Based on Sakai and Imiya 2009, "Fast Spectral Clustering with Random Projection and Sampling," algorithm 3.
      % Uses only landmarks as in the Nystrom approximation (see Platt 2005, "FastMap, MetricMap, and Landmark MDS are all Nystrom Algorithms"), but computes eigenvectors from number_landmarks whole columns instead of just the first number_landmarks rows of those columns.
      % NOTE
      % NOTE
      % NOTE
      % NOTE
      % NOTE
      % Should eventually implement this using Brand 2002's incremental SVD so all distances get used!
      
      % error('Not yet implemented!')
      % warning('Should respect options.desired_shape_space_dimensionality!')
      
      
      % if ~isempty(options.distances_to_compute)
        % error('method "%s" does not support option distances_to_compute being nonempty', options.method)
      % end
      
      % Obtain A and B (variable names from Platt 2005):
      landmark_distances = full(distances(landmark_indices, landmark_indices));
      remainder_distances = full(distances(~landmark_indices, landmark_indices));
      
      % Center A and B as in Landmark MDS:
      % A, Platt 2005, equation (13):
      landmark_centering_matrix = eye(number_landmarks) - 1 / number_landmarks;
      % This becomes landmark_centering_matrix * [E_ij^2 - 1/m * sum_q E_iq^2] = [(E_ij^2 - 1/m * sum_q E_iq^2) - sum_p (E_pj^2 - 1/m * sum_q E_pq^2)]?
      landmark_kernel = -1 / 2 * landmark_centering_matrix * (landmark_distances.^2) * landmark_centering_matrix;
      % B, Platt 2005, equation (15) (could change this to (14) at some point):
      % remainder_kernel = -1 / 2 * (remainder_distances.^2 - repmat(mean(remainder_distances.^2, 2), [1, number_landmarks]));
      remainder_kernel = -1 / 2 * (remainder_distances.^2 - repmat(mean(landmark_distances.^2, 1), [number_data - number_landmarks, 1]));
      
      % columns_kernel = [landmark_kernel; remainder_kernel];
      columns_kernel = zeros(number_data, size(landmark_kernel, 2));
      columns_kernel(landmark_indices, :) = landmark_kernel;
      columns_kernel(~landmark_indices, :) = remainder_kernel;

      % Don't need to do this because we already have a kernel matrix? distances isn't a similarity matrix... I assume I can do this because Sakai and Imiya 2009 notes that this is similar to the Nystrom approximation/extension.
      if false
        % Normalize affinity/weight matrix to create a kernel?
        % Sakai and Imiya 2009, equation (16), but with the square root of the pseudoinverse taken:
        row_sums = columns_kernel * ones(size(columns_kernel, 2), 1);
        row_sums(row_sums > 0) = 1 ./ sqrt(row_sums(row_sums > 0));
        diagonal_row_sums_pseudoinverse = spdiags(row_sums, 0, size(columns_kernel, 1), size(columns_kernel, 1));
        % Sakai and Imiya 2009, equation (17), but with the square root of the pseudoinverse taken:
        column_sums = ones(1, size(columns_kernel, 1)) * columns_kernel;
        column_sums(column_sums > 0) = 1 ./ sqrt(column_sums(column_sums > 0));
        size(column_sums), size(columns_kernel)
        keyboard
        diagonal_column_sums_pseudoinverse = spdiags(column_sums.', 0, size(columns_kernel, 2), size(columns_kernel, 2));
        % Sakai and Imiya 2009, equation (15):
        columns_kernel = diagonal_row_sums_pseudoinverse * columns_kernel * diagonal_column_sums_pseudoinverse;
      end
      
      [columns_kernel_left_singular_vectors, columns_kernel_singular_values, columns_kernel_right_singular_vectors] = svd(columns_kernel, 'econ');
      
      normalize_rows = false;
      % normalize_rows = true;
      normalize_rows_after_setting_dimensionality = false;
      % normalize_rows_after_setting_dimensionality = true;
      
      if normalize_rows_after_setting_dimensionality
        columns_kernel_left_singular_vectors = columns_kernel_left_singular_vectors(:, 1:min(end, desired_shape_space_dimensionality));
      end
      
      if normalize_rows
        % Sakai and Imiya 2009, algorithm 3, step 7:
        number_left_singular_vectors = size(columns_kernel_left_singular_vectors, 2)
        columns_kernel_left_singular_vectors = columns_kernel_left_singular_vectors ./ repmat(sqrt(sum(columns_kernel_left_singular_vectors.^2, 2)), [1, size(columns_kernel_left_singular_vectors, 2)]);
      end
      
      positions = columns_kernel_left_singular_vectors(:, 1:min(end, desired_shape_space_dimensionality));
      
      mass_matrix = speye(size(positions, 2));
      
      if false
        % Debug info:
        whos positions
        visible_state = get(gcf, 'Visible');
        set(gcf, 'Visible', 'on');
        scatter(positions(:, 1), positions(:, 2))
        axis image
        title(sprintf('In %s', mfilename), 'Interpreter', 'none')
        beep
        % keyboard
        pause
        set(gcf, 'Visible', visible_state);
      end

      % keyboard

      % error('Implementation yet unfinished!')
      
      
      if false
      
        if false
          % Debug info:
          reconstructed_landmark_kernel = landmark_kernel_eigenvectors * landmark_kernel_eigenvalues * landmark_kernel_eigenvectors'
          reconstructed_landmark_kernel_difference = reconstructed_landmark_kernel - landmark_kernel
          beep, keyboard
          % pause
        end
        
        % Sort by decreasing absolute eigenvalue:
        [~, landmark_kernel_eigenvalue_order] = sort(abs(diag(landmark_kernel_eigenvalues)), 1, 'descend');
        landmark_kernel_eigenvalues = landmark_kernel_eigenvalues(landmark_kernel_eigenvalue_order, landmark_kernel_eigenvalue_order);
        landmark_kernel_eigenvectors = landmark_kernel_eigenvectors(:, landmark_kernel_eigenvalue_order);
        
        % landmark_kernel_eigenvalues, pause
        
        % if options.use_schur_decomposition
          % if false
            % landmark_kernel_eigenvalues
            % warning('schur sometimes returns non-diagonal eigenvalues!')
          % end
          % landmark_kernel_eigenvalues = diag(diag(landmark_kernel_eigenvalues));
          % if false
            % landmark_kernel_eigenvalues
          % end
        % end
        
        
        % Remove unnecessary dimensions:
        
        
        unfiltered_landmark_kernel_eigenvalues = diag(landmark_kernel_eigenvalues).';
        embedding_info.unfiltered_landmark_kernel_eigenvalues = unfiltered_landmark_kernel_eigenvalues;
        
        good_eigenvalues = abs(unfiltered_landmark_kernel_eigenvalues) > eps(max(abs(unfiltered_landmark_kernel_eigenvalues))) * 1e2;
        
        if options.force_positive_definiteness
          % landmark_kernel_eigenvalues
          % unfiltered_landmark_kernel_eigenvalues
          unfiltered_eigenvalues_negative_proportion = abs(sum(unfiltered_landmark_kernel_eigenvalues(unfiltered_landmark_kernel_eigenvalues < 0))) / sum(abs(unfiltered_landmark_kernel_eigenvalues));
          
          % fprintf('Removing negative eigenvalues (%.3f of "variance")\n', unfiltered_eigenvalues_negative_proportion)
          good_eigenvalues = good_eigenvalues & (unfiltered_landmark_kernel_eigenvalues > 0);
          
          embedding_info.unfiltered_eigenvalues_negative_proportion = unfiltered_eigenvalues_negative_proportion;
        end
        
        % Respect desired_shape_space_dimensionality:
        good_eigenvalues(desired_shape_space_dimensionality + 1:end) = false;
        
        % Apply eigenvalue filtering:
        landmark_kernel_eigenvalues = diag(unfiltered_landmark_kernel_eigenvalues .* good_eigenvalues);
        
        eigenvalue_reordering = [find(good_eigenvalues), find(~good_eigenvalues)];
        landmark_kernel_eigenvalues = landmark_kernel_eigenvalues(eigenvalue_reordering, eigenvalue_reordering);
        landmark_kernel_eigenvectors = landmark_kernel_eigenvectors(:, eigenvalue_reordering);
        good_eigenvalues = good_eigenvalues(eigenvalue_reordering);
        landmark_kernel_eigenvalues = landmark_kernel_eigenvalues(good_eigenvalues, good_eigenvalues);
        landmark_kernel_eigenvectors = landmark_kernel_eigenvectors(:, good_eigenvalues);

        % S hat star, Schleif and Gisbrecht 2013, equation (4) (could use equation (5) at some point to work on distances directly?):
        columns_kernel = [landmark_kernel; remainder_kernel];
        % complete_kernel = columns_kernel * landmark_kernel_eigenvectors * pinv(landmark_kernel_eigenvalues) * landmark_kernel_eigenvectors.' * columns_kernel.';
        positions = columns_kernel * landmark_kernel_eigenvectors * sqrt(pinv(landmark_kernel_eigenvalues));
        
        mass_matrix = eye(size(positions, 2));
        
        if false
          % Debug info:
          if size(positions, 2) < 2
            beep
            keyboard
          end
        end
        
      end

    case 'optspace' %grj 10/30/14 - untested
        %can't force positive semi-definiteness
        
        tol = 1e-8;
        
        distances(isnan(distances)) = 0;
        
        [X S Y dist] = OptSpace(sparse(distances),[],[],tol);
        
        embedding_info.dist_complete = X*S*Y';
        mass_matrix = eye(size(positions,2));
        
    
    case 'svt_candes' %grj 10/30/14 - untested
        
        n = size(distances,1);
        Omega = find(~isnan(distances));
        b = distances(~isnan(distances));
        
        tau = 5*n; %default from Candes et al.
        delta = 2; %default
        
        
        [U, Sigma, V, numiter, out] = SVT(n, Omega, b, tau, delta);
        
        dist_complete = U*Sigma*V';
        
    case 'shortest path' %grj 11/2/14
        
        
        distances = sparse(distances);
        distances(isnan(distances)) = 0;
        dist = graphallshortestpaths(sparse(distances));
        
%         if embedding_info.weight_factor == 0;
%             positions = cmdscale(dist);
%         else
        [positions, stress] = mdscale(dist, min([size(dist,2), desired_shape_space_dimensionality]), 'Criterion', 'metricstress' , 'Weight', 1./(dist.^embedding_info.weight_factor), 'Replicates', embedding_info.replicates, 'Options', embedding_info.statset);
%         end
        
        embedding_info.stress = stress;
        embedding_info.positions_complete = positions;
        
        e = eig(positions*positions');
        e = sort(e, 'descend');
        e_pos = e.*(e>0);
        
        numdims = [find((cumsum(e_pos)/sum(e_pos)) >= embedding_info.explained_variance_threshold)];
        
        numdims = min(desired_shape_space_dimensionality, numdims);
        
        positions = positions(:,1:numdims(1));
        
        embedding_info.distances_complete = dist;
        embedding_info.eig = e;
        
        mass_matrix = eye(size(dist));
    case 'regress'

        a = distances;
        a_nans = isnan(a);
        
        weights = 1./(a.^embedding_info.weight_factor);
        weights(a_nans) = 0;
        
        a(a_nans) = inf;
        
        [positions, stress] = mdscale(a, min([size(a,2), desired_shape_space_dimensionality]), 'Criterion', 'metricstress', 'Weight', weights, 'Replicates', embedding_info.replicates, 'Options', embedding_info.statset, 'Start', embedding_info.Start);
    
        
        embedding_info.stress = stress;
        mass_matrix = eye(size(distances));
    case 'regress landmark'
        
        
        a = distances;
        %only use fully-populated columns
        a( any(isnan(a),1), any(isnan(a),1)) = nan;
        a_nans = isnan(a);
        
        a(logical(diag(ones(length(a),1)))) = 0;
        
        weights = 1./(a.^embedding_info.weight_factor);
%         weights(a_nans) = 0;
        
%         a(a_nans) = inf;
        
        [positions, stress] = mdscale(a,  min([size(a,2), desired_shape_space_dimensionality]), 'Criterion', 'metricstress', 'Weight', weights, 'Replicates', embedding_info.replicates, 'Options', embedding_info.statset, 'Start', embedding_info.Start);

        embedding_info.stress = stress;
        mass_matrix = eye(size(distances));
        
    case 'most_complete'
        %taken from complete_shape_space.m
        %grj 12/22/14
        
       %iteratively remove rows and columns with nans until we have a
        %nan-free matrix
        
        
        distances_incomplete_temp = distances;
        
        [ keepinds, distances_complete ] = get_complete_distance_matrix( distances_incomplete_temp );
        
       
         weights = 1./(distances_complete.^embedding_info.weight_factor);
%         if embedding_info.weight_factor == 0
%             [positions, e] = cmdscale( distances_complete );
%         else
            [positions, stress] = mdscale(distances_complete, min([size(distances_complete,2), desired_shape_space_dimensionality]), 'Weights', weights, 'Replicates', embedding_info.replicates, 'Options', embedding_info.statset);
            e = eig(positions*positions');
            e = e(end:-1:1);
%         end
        embedding_info.stress = stress;
        explained_variances = cumsum(e .* (e >= 0));
        explained_variances = explained_variances ./ explained_variances(end);
        sufficient_dimensionalities = find(explained_variances >= embedding_info.explained_variance_threshold); 

        numdims = min(sufficient_dimensionalities(1),desired_shape_space_dimensionality);

        positions = positions(:,1:numdims);

        positions_temp = nan(size(distances,1), size(positions,2));
        positions_temp(keepinds,:) = positions;
        
        positions = positions_temp;
        
        
        embedding_info.eig = e;
        
%     end
    
%     if numdims > size(positions, 2)
%       warning(['numdims = ', num2str(numdims), ' > size(positions, 2) = ', num2str(size(positions, 2))])
%       numdims = size(positions, 2);
%     end
    
    case 'mds'
         %taken from complete_shape_space.m
        %grj 12/22/14
        
       %iteratively remove rows and columns with nans until we have a
        %nan-free matrix
        
        
        distances_incomplete_temp = distances;
        
        [ keepinds, distances_complete ] = get_complete_distance_matrix( distances_incomplete_temp );
        
        [positions, e] = cmdscale( distances_complete );
        
        explained_variances = cumsum(e .* (e >= 0));
        explained_variances = explained_variances ./ explained_variances(end);
        sufficient_dimensionalities = find(explained_variances >= embedding_info.explained_variance_threshold); 

        numdims = min(sufficient_dimensionalities(1),desired_shape_space_dimensionality);

        positions = positions(:,1:numdims);

        positions_temp = nan(size(distances,1), size(positions,2));
        positions_temp(keepinds,:) = positions;
        
        positions = positions_temp;
        
        
        embedding_info.eig = e;
  
          
         
    otherwise
      error('Unrecognized method "%s"', embedding_info.method)
  end

  
  if false
    % Debug info:
    visible_state = get(gcf, 'Visible');
    set(gcf, 'Visible', 'on');
    scatter(positions(:, 1), positions(:, 2))
    axis image
    title(sprintf('In %s', mfilename), 'Interpreter', 'none')
    beep
    % keyboard
    pause
    set(gcf, 'Visible', visible_state);
  end
  
  
  positions_temp = nan(length(keepinds_master), size(positions,2));
  positions_temp(keepinds_master,:) = positions;
  
  positions = positions_temp;
  
  embedding_info.total_wall_time = toc(embedding_info.total_wall_time);
  embedding_info.total_cpu_time = cputime - embedding_info.total_cpu_time;
  
  
end

