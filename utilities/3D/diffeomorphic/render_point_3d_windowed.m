function r_structure = render_point_3d_windowed(...
  point, shape_space, image_function, epsilon, ...
  options)
  % Synthesize a shape for a point in shape_space. image_function returns the shape image corresponding to a point. Interpolates between shape images from the same Delaunay cell as point.
  %
  % 2013-02-11 tebuck: if shape_space is a cell array of length 4, use the fourth argument (the distances matrix used to compute the shape space) to do faster interpolations.
  % 03/30/2016 xruan: set options for using faster code
  
  options.just_compute_distance = false;
  % xruan 03/30/2016
  options = ml_initparam(options, struct( ...
    'distance_computing_method', 'faster' ...
    ));
  distance_computing_method = options.distance_computing_method;
  
  [M, N, P] = size(image_function(1));
  window_size = options.window_radius * [1, 1, 0] + [0, 0, P];

  
  Y2 = shape_space{1};
  convex_hull = shape_space{2};
  tes = shape_space{3};

  use_distance_matrix = false;
  if length(shape_space) > 3
    use_distance_matrix = true;
    distance_matrix = shape_space{4};
  end

  %Locate sample point in triangle mesh
  %Y2, tes, point
  simplex_idx = tsearchn(Y2,tes,point);
  xidx = tes(simplex_idx,:);
  xtri = Y2(xidx,:);

  distances = sqrt( sum( (xtri-repmat( point, [length(xidx),1])).^2, 2 ));
  [minv,idx] = min( distances );

  compression_method = 'none';
  if isfield(options, 'use_compression') && options.use_compression
    compression_method = 'gzip';
  end
  
  if distances(idx)<= epsilon
    % r_structure.interpolated_image = window_compress_image(...
      % image_function(xidx(idx)), ...
      % window_size);
    r_structure.interpolated_image = WindowedImage(...
      image_function(xidx(idx)), window_size, compression_method);
    r_structure.total_wall_time = 0;
    r_structure.total_cpu_time = 0;
    return
  end

  % Shape interpolation path
  [lambda, simplex_vertex_ordering] = shape_interp_track(xtri',point');
  simplex_vertex_ordering_temp = simplex_vertex_ordering; %save this because elsewhere it may change
  
  shape_ordering = xidx(simplex_vertex_ordering);
  
  
  
  % % Get distances corresponding to each lambdas:
  if use_distance_matrix
    simplex_distance_matrix = distance_matrix(shape_ordering, shape_ordering)
    % whos lambda ordering
    % % Find (approximate) true distances of interpolated points by using MDS on the simplex's vertices, interpolating between them by given amounts, and recomputing the lower-dimensional simplex's distance matrix from that:
    % simplex_distance_matrix = distance_matrix(ordering, ordering)
    % lambdas_as_distances = zeros(1, length(ordering));
    % % keyboard
    % for distance_index = 1:length(lambdas_as_distances)
      % [simplex_mds, e] = cmdscale(simplex_distance_matrix);
      % simplex_distance_matrix
      % lambdas_as_distances(distance_index) = distance_matrix()
    % end
  end
  % lambda, ordering

  % c_0 = x_1 from the paper:
  temp = image_function(shape_ordering(end));
  total_wall_time = 0;
  total_cpu_time = 0;
  
  ri = cell(1, length(lambda));
  
  indorder = length(lambda):-1:1;
  
  for i=indorder
    % lambda_1 (the last one in terms of array index) = (x_2 - c_1) / (c_1 - x_1),
    % (x_2 - c_1) = lambda_1 (c_1 - x_1),
    % so the distance between x_1 and c_1 is 

    % lambda_k:
    amount = lambda(i);
    % i
    % shape_ordering(i)
    % x_(k+1):
    start = image_function(shape_ordering(i));
    % c_(k-1):
    last = temp;
    if isfield(options, 'registration_error_function')
      options.convergence_registration_error = ...
          options.registration_error_function(start, last);
    end

    if use_distance_matrix
      % whos lambda simplex_vertex_ordering shape_ordering
      % Find (approximate) true distances of interpolated points (to speed up image interpolation) by using MDS on the simplex's vertices, interpolating between them by given amounts, and recomputing the lower-dimensional simplex's distance matrix from that:
      % simplex_distance_matrix = distance_matrix(shape_ordering, shape_ordering)
      % [simplex_mds, e] = cmdscale(simplex_distance_matrix)
      [simplex_mds, e] = mdscale(simplex_distance_matrix, size(simplex_distance_matrix, 1) - 1)
      simplex_mds_distances = squareform(pdist(simplex_mds))
      vertex_indices = [i, length(simplex_vertex_ordering)]
      % Get the indices of the two shapes being interpolated:
      temp_vertex_ordering = simplex_vertex_ordering(vertex_indices);
      temp_shape_ordering = shape_ordering(vertex_indices);
      % Get these shapes' distances:
      % temp_vertices = simplex_mds(temp_vertex_ordering, :);
      temp_vertices = simplex_mds(vertex_indices, :);
      % Get these shapes' (approximate) true distances:
      % temp_distances = simplex_distance_matrix(temp_vertex_ordering, :);
      temp_distances = simplex_distance_matrix(vertex_indices, :);
      % Clear those indices from the non-temporary variables as they will be replaced by values related to the interpolated point:
      % simplex_distance_matrix(temp_vertex_ordering, temp_vertex_ordering) = [];
      % simplex_distance_matrix(temp_vertex_ordering, :) = [];
      % simplex_distance_matrix(:, temp_vertex_ordering) = [];
      simplex_distance_matrix(vertex_indices, :) = [];
      simplex_distance_matrix(:, vertex_indices) = [];
      simplex_vertex_ordering(vertex_indices) = [];
      shape_ordering(vertex_indices) = [];
      % simplex_mds(temp_vertex_ordering, :) = [];
      simplex_mds(vertex_indices, :) = [];
      % amount, temp_vertex_ordering, temp_shape_ordering
      % keyboard
      % Create the interpolated vertex:
      interpolated_vertex = temp_vertices(1, :) * amount / (1 + amount) + temp_vertices(2, :) * 1 / (1 + amount);
      simplex_mds(end + 1, :) = interpolated_vertex;
      simplex_vertex_ordering(end+1) = length(simplex_vertex_ordering) + 1;
      simplex_distance_matrix(end+1, 1:end) = 0;
      simplex_distance_matrix(1:end, end+1) = 0;
      simplex_distance_matrix(1:end, end) = pdist2(simplex_mds, interpolated_vertex);
      simplex_distance_matrix(end, 1:end) = simplex_distance_matrix(1:end, end)';
      % This value is never used:
      shape_ordering(end+1) = 0;
      % 
      % known_distance = temp_distances(1, temp_vertex_ordering(2));
      known_distance = temp_distances(1, vertex_indices(2));
      stop_source_distance = known_distance * 1 / (1 + amount);
      known_distance, stop_source_distance
      options.known_distance = known_distance;
      options.stop_source_distance = stop_source_distance;
      options.lambdas_are_distances_from_source = true;
      diary_state = get(0, 'Diary'); diary('off'); diary(diary_state);
      % keyboard
      % fprintf('>>>> HACK\n')
      % xruan 03/30/2016
      if strcmp(distance_computing_method, 'faster')
          ri{i} = Greedy3D_lambda_pre(...
            start, last, stop_source_distance, options);
          temp = ri{i}.interpolated_image;
      else
          ri{i} = Greedy3D_lambda_pre_compressed(...
            start, last, stop_source_distance, options);
          temp = ri{i}.interpolated_image.get_data();
      end
    else
      if strcmp(distance_computing_method, 'faster')
          ri{i} = Greedy3D_lambda_pre(...
            start, last, amount, options);
          temp = ri{i}.interpolated_image;
      else
          ri{i} = Greedy3D_lambda_pre_compressed(...
            start, last, amount, options);
        temp = ri{i}.interpolated_image.get_data();
      end        
    end
    % fprintf('>>>> HACK\n')
    
    total_wall_time = total_wall_time + ri{i}.total_wall_time;
    total_cpu_time = total_cpu_time + ri{i}.total_cpu_time;
  end

  %image = temp;
  % r_structure.interpolated_image = window_compress_image(...
    % temp, window_size);
  r_structure.ri = ri(indorder);
  r_structure.lambda = lambda(indorder);
  r_structure.point = point;
  r_structure.simplex_vertex_ordering = simplex_vertex_ordering_temp(end:-1:1);
  r_structure.xidx = xidx;
  r_structure.xtri = xtri;
  
  r_structure.interpolated_image = WindowedImage(...
    temp, window_size, compression_method);
  r_structure.total_wall_time = total_wall_time;
  r_structure.total_cpu_time = total_cpu_time;
