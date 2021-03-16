function [result_distances] = mass_distance_function(given_vectors1, given_vectors2, given_mass_matrix)
  % If first argument gives multiple vectors, distances are between corresponding rows, not a distance matrix:
  if ~exist('given_mass_matrix', 'var')
    given_mass_matrix = eye(size(given_vectors1, 2));
  end
  if size(given_vectors1, 1) == 1
    given_vectors1 = repmat(given_vectors1, [size(given_vectors2, 1), 1]);
  end
  differences = given_vectors1 - given_vectors2;
  try
    result_distances = sqrt(sum(differences .* (given_mass_matrix * differences.').', 2));
  catch
    whos given_vector* given_mass_matrix result_distances
    keyboard
    error
  end
end
