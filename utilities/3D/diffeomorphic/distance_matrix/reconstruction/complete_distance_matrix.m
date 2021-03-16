function [distances] = complete_distance_matrix(distances, options)
  % Provides a consistent interface to complete a partially observed distance matrix using one of several methods. The observations must have a certain structure for certain methods, and those should be documented here...
  %
  % Methods nystrom-euclidean, nystrom-dissimilarity-clip, nystrom-dissimilarity-orthogonalize: Must observe the first desired_shape_space_dimensionality + 2 columns.
  %   
  % 2013-11-12 tebuck: copied from register_shapes_to_mean_shape.m.
  % % 2013-11-24 tebuck: adding nystrom approximation-inspired, missing-data tolerant, linear-time .
  

  default_options = struct();
  default_options.desired_shape_space_dimensionality = inf;
  % If this function is given and the method correctly selected, new distances can be computed upon request (this function should take care of things like using chunk_start, and embed_partial_distance_matrix should use its own random number generator instance instead of the default one):
  % default_options.distance_function = @(given_index)[];
  default_options.distance_function = [];
  default_options.method = 'direct';

  if ~exist('options', 'var')
    options = default_options; 
  else
    options = process_options_structure(default_options, options);
  end


  % Copy some options to local variables:
  desired_shape_space_dimensionality = options.desired_shape_space_dimensionality;

  
  number_data = size(distances, 1);
  
  
  
  % Find landmarks using already observed entries:
  
  use_fixed_dimensionality_landmarks = false;
  % use_fixed_dimensionality_landmarks = true;
  if use_fixed_dimensionality_landmarks
    % The number of landmarks is chosen based on a Euclidean matrix's rank if its points have dimensionality desired_shape_space_dimensionality (see Drineas et al. 2006, "Distance Matrix Reconstruction from Incomplete Distance Information for Sensor Network Localization"):
    number_landmarks = desired_shape_space_dimensionality + 2;
    landmark_indices = [true(1, number_landmarks), false(1, number_data - number_landmarks)];
  else
    % Number of observations not on the diagonal (assumes no duplicate points):
    column_number_observations = sum(~isnan(distances) & distances ~= 0);
    landmark_indices = column_number_observations == max(column_number_observations);
    number_landmarks = sum(landmark_indices);
  end
  % keyboard
  
  if number_landmarks == 0
    error('No landmarks observed! No methods currently complete randomly sampled distance matrices.')
  end
  
  
  
  % Reconstruct entries with nans:
  switch options.method
    case 'direct'
    
    
      distances = reconstruct_distance_matrix(distances, desired_shape_space_dimensionality);
      
      
    case 'mishra'
    
    
      blanked_entries = isnan(distances);
      % [rows, columns] = find(~squareform(blanked_entries));
      [rows, columns] = find(~blanked_entries);
      good_indices = rows > columns;
      rows = rows(good_indices);
      columns = columns(good_indices);
      % given_distances = distances(~blanked_entries);
      given_distances = distances(good_indices);
      mishra_parameters = struct();
      mishra_parameters.rmax = desired_shape_space_dimensionality + 2;
      mishra_parameters.verb = false;
      % whos rows columns given_distances
      % number_data, desired_shape_space_dimensionality
      % [mishra_data, info_structure] = lowrank_dist_completion(@tr_dist_completion, rows, columns, given_distances', randn(number_data, desired_shape_space_dimensionality + 2), mishra_parameters);
      [mishra_data, info_structure] = lowrank_dist_completion(@tr_dist_completion, rows, columns, given_distances, randn(number_data, desired_shape_space_dimensionality + 2), mishra_parameters);
      warning('Why was pdist output originally squared???')
      % distances = pdist(mishra_data)'.^2;
      distances = pdist(mishra_data)';
      distances = squareform(distances);
      
      
    case 'nystrom-euclidean'
    
    
      % Use only landmarks in the Nystrom approximation (see Platt 2005, "FastMap, MetricMap, and Landmark MDS are all Nystrom Algorithms"):
      
      % if ~isempty(options.distances_to_compute)
        % error('method "%s" does not support option distances_to_compute being nonempty', options.method)
      % end
      
      % Obtain A and B (variable names from Platt 2005):
      landmark_distances = distances(landmark_indices, landmark_indices);
      remainder_distances = distances(~landmark_indices, landmark_indices);
      % % Use the orientation in the paper:
      % remainder_distances = distances(landmark_indices, ~landmark_indices);
      
      % Center A and B as in Landmark MDS:
      % A, Platt 2005, equation (13):
      landmark_centering_matrix = eye(number_landmarks) - 1 / number_landmarks;
      % This becomes landmark_centering_matrix * [E_ij^2 - 1/m * sum_q E_iq^2] = [(E_ij^2 - 1/m * sum_q E_iq^2) - sum_p (E_pj^2 - 1/m * sum_q E_pq^2)]?
      landmark_kernel = -1 / 2 * landmark_centering_matrix * (landmark_distances.^2) * landmark_centering_matrix;
      % B, Platt 2005, equation (13):
      remainder_kernel = -1 / 2 * (remainder_distances.^2 - repmat(mean(landmark_distances.^2, 2).', [number_data - number_landmarks, 1]));
      
      
      if false
        % Debug info:
        % Check the matrix multiplication version with this one:
        landmark_squared_distance_row_means = mean(landmark_distances.^2, 2);
        landmark_squared_distance_column_means = mean(landmark_distances.^2, 1);
        landmark_squared_distance_mean = mean(landmark_distances(:).^2);
        for i = find(landmark_indices)
          for j = find(landmark_indices)
            landmark_kernel2(i, j) = -1 / 2 * (landmark_distances(i, j).^2 - landmark_squared_distance_column_means(j) - landmark_squared_distance_row_means(i) + landmark_squared_distance_mean);
          end
        end
        landmark_kernel_error = norm(landmark_kernel - landmark_kernel2)
        
        for i = find(~landmark_indices)
          for j = find(landmark_indices)
            remainder_kernel2(i, j) = -1 / 2 * (remainder_distances(i, j).^2 - landmark_squared_distance_row_means(j));
          end
        end
        remainder_kernel_error = norm(remainder_kernel - remainder_kernel2)
      end
      
      
      % % Don't do the MDS part of this here:
      % % Decompose A, Platt 2005, equation (7):
      % [landmark_kernel_eigenvectors, landmark_kernel_eigenvalues] = schur(landmark_kernel);
      
      
      % if false
        % % Debug info:
        % % Check the approximation:
        % decomposition_error = norm(U * Gamma * U' - landmark_kernel)
      % end
      
      
      % K tilde, Platt 2005, equation (11):
      complete_kernel = [landmark_kernel, remainder_kernel.'; remainder_kernel, remainder_kernel * pinv(landmark_kernel) * remainder_kernel.'];
      % Make symmetric:
      complete_kernel = 1 / 2 * (complete_kernel + complete_kernel.');
      % For Euclidean distances of low enough dimensionality, this should work, otherwise there might be complex values:
      warning('Should this use equation (3) in Belongie et al. 2002 for efficiency?')
      complete_distances = sqrt(repmat(diag(complete_kernel), [1, number_data]) + repmat(diag(complete_kernel), [1, number_data]).' - 2 * complete_kernel);
      warning('>>>> HACK, setting imaginary distances to zero'), complete_distances = real(complete_distances);
      distances = complete_distances;
      
      
    case 'nystrom-dissimilarity-clip'
    
    
      % Based on method 'nystrom-euclidean' but using eigenvalue correction from Schleif and Gisbrecht 2013, "Data analysis of (non-)metric proximities at linear costs:"
      
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
      % B, Platt 2005, equation (13):
      remainder_kernel = -1 / 2 * (remainder_distances.^2 - repmat(mean(landmark_distances.^2, 2).', [number_data - number_landmarks, 1]));
      
      
      % Decompose A, Platt 2005, equation (7):
      [landmark_kernel_eigenvectors, landmark_kernel_eigenvalues] = schur(landmark_kernel);
      
      % Correct eigenvalues:
      % Clip:
      landmark_kernel_eigenvalues(landmark_kernel_eigenvalues < 0) = 0;
      
      
      % S hat star, Schleif and Gisbrecht 2013, equation (4) (could use equation (5) at some point to work on distances directly?):
      columns_kernel = [landmark_kernel; remainder_kernel];
      complete_kernel = columns_kernel * landmark_kernel_eigenvectors * pinv(landmark_kernel_eigenvalues) * landmark_kernel_eigenvectors.' * columns_kernel.';
      % Make symmetric:
      complete_kernel = 1 / 2 * (complete_kernel + complete_kernel.');

      % This should work in general according to Schleif and Gisbrecht 2013 (right before (5)):
      complete_distances = sqrt(repmat(diag(complete_kernel), [1, number_data]) + repmat(diag(complete_kernel), [1, number_data]).' - 2 * complete_kernel);
      % warning('>>>> HACK, setting imaginary distances to zero'), complete_distances = real(complete_distances);
      distances = complete_distances;
      
      
      % Should instead get distances from Euclidean embedding?
      
      
      % keyboard
      % error
      
      
    case 'nystrom-dissimilarity-orthogonalize'
    
    
      % Based on method 'nystrom-euclidean' but using orthogonalization from Belongie et al. 2002, "Spectral Partitioning with Indefinite Kernels Using the Nystrom Extension:"
      
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
      % B, Platt 2005, equation (13):
      remainder_kernel = -1 / 2 * (remainder_distances.^2 - repmat(mean(landmark_distances.^2, 2).', [number_data - number_landmarks, 1]));
      
      
      % Decompose A, Platt 2005, equation (7):
      [landmark_kernel_eigenvectors, landmark_kernel_eigenvalues] = schur(landmark_kernel);
      

      % Correct eigenvalues:
      % should_clip_eigenvalues = false
      should_clip_eigenvalues = true
      if should_clip_eigenvalues
        % Clip:
        minimum_landmark_kernel_eigenvalue = min(diag(landmark_kernel_eigenvalues))
        large_negative_eigenvalue_threshold = -mean(abs(diag(landmark_kernel_eigenvalues))) * eps * 1e2
        if any(diag(landmark_kernel_eigenvalues) < large_negative_eigenvalue_threshold)
          warning('Clipping large negative eigenvalues!')
        end
        landmark_kernel_eigenvalues(landmark_kernel_eigenvalues < 0) = 0;
      end
      
      
      % Final formula for W hat, Belongie et al. 2002, section 4:
      % columns_kernel = [landmark_kernel; remainder_kernel];
      % complete_kernel = columns_kernel * landmark_kernel_eigenvectors * pinv(landmark_kernel_eigenvalues) * landmark_kernel_eigenvectors.' * columns_kernel.';
      % columns_kernel = [landmark_kernel; remainder_kernel * landmark_kernel_eigenvectors * pinv(landmark_kernel_eigenvalues)];
      projected_columns_kernel = [landmark_kernel_eigenvectors; remainder_kernel * landmark_kernel_eigenvectors * pinv(landmark_kernel_eigenvalues)];
      positions = projected_columns_kernel * sqrt(landmark_kernel_eigenvalues);
      [positions_eigenvectors, positions_eigenvalues] = schur(positions.' * positions);
      complete_kernel_eigenvectors = positions * positions_eigenvectors * sqrt(pinv(positions_eigenvalues));
      complete_kernel = complete_kernel_eigenvectors * positions_eigenvalues * complete_kernel_eigenvectors.';
      % Make symmetric:
      complete_kernel = 1 / 2 * (complete_kernel + complete_kernel.');

      
      if false
        % Debug info:
        % Check the above for its requirements:
        whos columns* positions* complete_kernel*
      end
      

      % This should work in general according to Schleif and Gisbrecht 2013 (right before (5)):
      complete_distances = sqrt(repmat(diag(complete_kernel), [1, number_data]) + repmat(diag(complete_kernel), [1, number_data]).' - 2 * complete_kernel);
      % warning('>>>> HACK, setting imaginary distances to zero'), complete_distances = real(complete_distances);
      distances = complete_distances;
      
      
      % Should instead get distances from Euclidean embedding?
      
      
      % keyboard
      % error
      
      
    otherwise
      error('Unrecognized method "%s"', options.method)
  end
  
  
end

