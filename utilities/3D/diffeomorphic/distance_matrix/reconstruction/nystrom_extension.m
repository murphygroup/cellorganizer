function [U, S] = nystrom_extension(kernel, options)
  % Provides a consistent interface to compute Nystrom extensions from a partially observed real symmetric matrix given_values using one of several methods. Some number of (the same) entire rows and columns in given_values should be provided.
  %
  %
  % Tests
  % =====
  % n = 100, m = ceil(sqrt(n) / 2), a = randn(n, m); a = a * a'; ac = spalloc(n, n, n * m * 2 - m^2); ac(:, 1:ceil(m * 1.5)) = a(:, 1:ceil(m * 1.5)); ac(1:ceil(m * 1.5), :) = a(1:ceil(m * 1.5), :); [au, as] = nystrom_extension(ac, struct('method', 'classical', 'force_positive_definiteness', false)); norm(a - au * as * au', 'fro') / norm(a, 'fro')
  %
  %
  % 2014-03-16 tebuck: Copied from embed_partial_distance_matrix.m.
  

  default_options = struct();
  default_options.method = 'classical';
  % Clip negative eigenvalues and such so the results have a positive definite distance matrix regardless of the input's definiteness:
  default_options.force_positive_definiteness = false;
  default_options.maximum_rank = inf;

  
  if ~exist('options', 'var')
    options = default_options; 
  else
    options = process_options_structure(default_options, options);
  end


  % Copy some options to local variables:
  maximum_rank = options.maximum_rank;

  
  number_data = size(kernel, 1);
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Common code:
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  given_entries = kernel ~= 0;
  given_entries_with_diagonals = given_entries | logical(speye(number_data));
  given_entries_proportion = full(sum(sum(given_entries))) ./ numel(given_entries);
  
  
  % Find landmarks using already observed entries:
  
  % Number of observations not on the diagonal (assumes no duplicate points):
  % column_number_observations = sum(~isnan(kernel) & kernel ~= 0);
  column_number_observations = sum(given_entries);
  landmark_indices = column_number_observations == max(column_number_observations);
  number_landmarks = sum(landmark_indices);

  
  if number_landmarks == 0
    error('No landmarks observed!')
  end
  
  
  
  function [result_eigenvectors, result_eigenvalues] = spectral_decomposition(given_matrix, given_sort_by_absolute_eigenvalues)
    if ~exist('given_sort_by_absolute_eigenvalues', 'var') || isempty(given_sort_by_absolute_eigenvalues)
      given_sort_by_absolute_eigenvalues = false;
    end
    [result_eigenvectors, result_eigenvalues] = eig(given_matrix);
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
  
  
  % Obtain A and B (variable names from Platt 2005):
  % Obtain A_M and B_M (variable names from Nemtsov):
  landmark_kernel = full(kernel(landmark_indices, landmark_indices));
  remainder_kernel = full(kernel(landmark_indices, ~landmark_indices));
  data_order_restoration = zeros(number_data, 1);
  data_order_restoration(landmark_indices) = (1:number_landmarks)';
  data_order_restoration(~landmark_indices) = (number_landmarks + 1:number_data)';
  
  % Decompose A_M into U and Lambda, Nemtsov et al. 2013, equation (2.1):
  [landmark_kernel_eigenvectors, landmark_kernel_eigenvalues] = spectral_decomposition(landmark_kernel, ~options.force_positive_definiteness);
  
  % Reconstruct entries with nans:
  switch options.method
    case 'classical'
      
      
      % error('Not yet implemented!')
      
      % Based on Nemtsov et al. 2013, "Matrix Compression using the Nystrom Method," section 2.1:
      
      % Get extended left eigenvectors Uhat from U and Utilde, Nemtsov et al. 2013, equation (2.4), (2.6):
      % kernel_left_eigenvectors = [landmark_kernel_eigenvectors; remainder_kernel' * landmark_kernel_eigenvectors * pinv(landmark_kernel_eigenvalues)];
      % kernel_left_eigenvectors = [landmark_kernel_eigenvectors; remainder_kernel' * landmark_kernel_eigenvectors * sparse(pinv(full(landmark_kernel_eigenvalues)))];
      landmark_kernel_eigenvalues_inverse = landmark_kernel_eigenvalues;
      landmark_kernel_eigenvalues_inverse(landmark_kernel_eigenvalues_inverse ~= 0) = 1 ./ landmark_kernel_eigenvalues_inverse(landmark_kernel_eigenvalues_inverse ~= 0);
      kernel_left_eigenvectors = [landmark_kernel_eigenvectors; remainder_kernel' * landmark_kernel_eigenvectors * landmark_kernel_eigenvalues_inverse];
      kernel_left_eigenvectors = kernel_left_eigenvectors(data_order_restoration, :);
      
      U = kernel_left_eigenvectors;
      S = landmark_kernel_eigenvalues;

      
    case 'belongie2002'
    
    
      error('Not yet implemented!')
    
    
      % Based on method 'classical' but using orthogonalization from Belongie et al. 2002, "Spectral Partitioning with Indefinite Kernels Using the Nystrom Extension:"
      
      
      
      
    case 'sdm-nystrom-kernel-orthogonalize'
    
    
      % error('Not yet implemented!')
      error('Not yet finished!')
    
      % Based on method 'nystrom-dissimilarity-orthogonalize' but uses the Nystrom extension on the squared distance matrix, converts that to a kernel by double centering, and orthogonalizes the kernel as in Belongie et al. 2002, "Spectral Partitioning with Indefinite Kernels Using the Nystrom Extension." Thus this is my novel approach inspired by Belongie et al. 2002 and Liu 2009, "Spectral mesh segmentation" (5.4):
      
      % Variable names from "Taraz Shape space estimation 004.docx":
      
      
      
    case 'approximate-evd'
    
    
      % error('Not yet implemented!')
      % error('Not yet finished!')
      
      
      % Based on Nemtsov et al. 2013, "Matrix Compression using the Nystrom Method," section 2.2.3:
      
      % Get extended left eigenvectors Uhat from U and Utilde, Nemtsov et al. 2013, equation (2.4), (2.6), (2.9):
      landmark_kernel_square_root_inverse = pinv(general_square_root(landmark_kernel));
      % Form the matrix G in section 2.2.3:
      kernel_scaled_matrix = [landmark_kernel; remainder_kernel'] * landmark_kernel_square_root_inverse;
      % Compute the scatter matrix's spectral decomposition U_S, Lambda_S:
      % [kernel_scaled_matrix_scatter_eigenvectors, kernel_scaled_matrix_scatter_eigenvalues] = spectral_decomposition(kernel_scaled_matrix' * kernel_scaled_matrix, ~options.force_positive_definiteness);
      [kernel_scaled_matrix_scatter_eigenvectors, kernel_scaled_matrix_scatter_eigenvalues] = spectral_decomposition(landmark_kernel + landmark_kernel_square_root_inverse * remainder_kernel * remainder_kernel' * landmark_kernel_square_root_inverse, ~options.force_positive_definiteness);
      % U_o:
      orthogonal_eigenvectors = kernel_scaled_matrix * kernel_scaled_matrix_scatter_eigenvectors * pinv(general_square_root(kernel_scaled_matrix_scatter_eigenvalues));
      orthogonal_eigenvalues = kernel_scaled_matrix_scatter_eigenvalues;
      
      orthogonal_eigenvectors = orthogonal_eigenvectors(data_order_restoration, :);
      
      if false
        fprintf('landmark_kernel reconstruction relative Frobenius error: %e\n', norm(landmark_kernel - kernel_corrected_left_eigenvectors(landmark_indices, :) * kernel_corrected_eigenvalues * kernel_corrected_left_eigenvectors(landmark_indices, :)', 'fro') ./ norm(landmark_kernel, 'fro'))
        fprintf('kernel_corrected_left_eigenvectors orthogonality relative Frobenius error: %e\n', norm(eye(size(kernel_corrected_left_eigenvectors, 2)) - kernel_corrected_left_eigenvectors' * kernel_corrected_left_eigenvectors, 'fro') ./ norm(eye(size(kernel_corrected_left_eigenvectors, 2)), 'fro'))
        fprintf('remainder_kernel reconstruction relative Frobenius error: %e\n', norm(remainder_kernel - kernel_corrected_left_eigenvectors(landmark_indices, :) * kernel_corrected_eigenvalues * kernel_corrected_left_eigenvectors(~landmark_indices, :)', 'fro') ./ norm(remainder_kernel, 'fro'))
        kernel_corrected_left_eigenvectors
        kernel_corrected_eigenvalues
        beep, keyboard
        % pause
      end
      
      
      U = orthogonal_eigenvectors;
      S = orthogonal_eigenvalues;
      
      
    case 'pekalska-pseudo-euclidean'
    
      
      error
      
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
      
      if false
        warning('>>>> HACK, not removing any eigenvalues!')
        good_eigenvalues = true(size(unfiltered_landmark_kernel_eigenvalues));
      else
        % % Zero eigenvalues:
        % good_eigenvalues = diag(landmark_kernel_eigenvalues) ~= 0;
        % Small absolute value eigenvalues:
        % good_eigenvalues = abs(diag(landmark_kernel_eigenvalues)) > eps(max(abs(diag(landmark_kernel_eigenvalues)))) * 1e2;
        good_eigenvalues = abs(unfiltered_landmark_kernel_eigenvalues) > eps(max(abs(unfiltered_landmark_kernel_eigenvalues))) * 1e2;
        
        if options.force_positive_definiteness
          % landmark_kernel_eigenvalues
          % unfiltered_landmark_kernel_eigenvalues
          unfiltered_eigenvalues_negative_proportion = abs(sum(unfiltered_landmark_kernel_eigenvalues(unfiltered_landmark_kernel_eigenvalues < 0))) / sum(abs(unfiltered_landmark_kernel_eigenvalues));
          
          % fprintf('Removing negative eigenvalues (%.3f of "variance")\n', unfiltered_eigenvalues_negative_proportion)
          good_eigenvalues = good_eigenvalues & (unfiltered_landmark_kernel_eigenvalues > 0);
        end
      end
      
      % Respect maximum_rank:
      good_eigenvalues(maximum_rank + 1:end) = false;
      
      % Apply eigenvalue filtering:
      landmark_kernel_eigenvalues = diag(unfiltered_landmark_kernel_eigenvalues .* good_eigenvalues);
      
      eigenvalue_reordering = [find(good_eigenvalues), find(~good_eigenvalues)];
      landmark_kernel_eigenvalues = landmark_kernel_eigenvalues(eigenvalue_reordering, eigenvalue_reordering);
      landmark_kernel_eigenvectors = landmark_kernel_eigenvectors(:, eigenvalue_reordering);
      good_eigenvalues = good_eigenvalues(eigenvalue_reordering);
      landmark_kernel_eigenvalues = landmark_kernel_eigenvalues(good_eigenvalues, good_eigenvalues);
      landmark_kernel_eigenvectors = landmark_kernel_eigenvectors(:, good_eigenvalues);
      
      
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
        if options.force_positive_definiteness
          euclidean_dimensions = sign(diag(mass_matrix)) > 0;
          positions = positions(:, euclidean_dimensions);
          mass_matrix = eye(size(positions, 2));
        end
        dimensionality_after_forcing = size(positions, 2);
      % end

      
      % keyboard
      % error
    
      
    otherwise
      error('Unrecognized method "%s"', options.method)
  end


  % Filter eigenvalues:
  
  % Sort by decreasing absolute eigenvalue:
  [~, S_order] = sort(abs(diag(S)), 1, 'descend');
  S = S(S_order, S_order);
  U = U(:, S_order);
  
  % Remove unnecessary dimensions:
  
  unfiltered_S = diag(S).';
  
  S_good = abs(unfiltered_S) > eps(max(abs(unfiltered_S))) * 1e2;
  
  if options.force_positive_definiteness
    % S
    % unfiltered_S
    unfiltered_eigenvalues_negative_proportion = abs(sum(unfiltered_S(unfiltered_S < 0))) / sum(abs(unfiltered_S));
    
    % fprintf('Removing negative eigenvalues (%.3f of "variance")\n', unfiltered_eigenvalues_negative_proportion)
    S_good = S_good & (unfiltered_S > 0);
  end
  
  % Respect maximum_rank:
  % S_good(maximum_rank + 1:end) = false;
  S_good(cumsum(S_good) > maximum_rank) = false;
  
  % Apply eigenvalue filtering:
  S = diag(unfiltered_S .* S_good);
  
  eigenvalue_reordering = [find(S_good), find(~S_good)];
  S = S(eigenvalue_reordering, eigenvalue_reordering);
  U = U(:, eigenvalue_reordering);
  S_good = S_good(eigenvalue_reordering);
  S = S(S_good, S_good);
  U = U(:, S_good);


  
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
  
  
  
end

