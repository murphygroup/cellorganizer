function [plane_parameters] = fit_plane_to_synapse_in_mesh(single_cell_mesh, synapse_point, options)
% function [plane_parameters] = fit_plane_to_synapse_in_mesh(single_cell_mesh, synapse_point, options)
% Dependencies:
%% GrTheory, toolbox_graph, toolbox_wavelet_meshes (all from File Exchange)
% toolbox_graph, toolbox_fast_marching, inpolyhedron, lowner (all from File Exchange)
% 2012-10-18 tebuck: Copied from segment_mesh_using_fast_marching.m.
% 2013-03-03 tebuck: Adding option to switch methods (synapse_detection_method). Added synapse_detection_method = 'synapse direction'. Fixed normal_smoothed by normalizing. Fixed neighbors_to_add by using acos instead of cos. Changed neighbors_to_add to use the best new neighbor instead of all within the threshold.
% 2013-03-27 tebuck: Completing synapse_detection_method = 'synapse direction'.

  default_options.curvature_smoothing = 2;
  default_options.debug_plot = false;
  % default_options.weighting_function = @(residuals, distances)1 ./ (distances + min(distances(distances > 0))).^2;
  default_options.weighting_function = @(residuals, distances)1 ./ (distances / mean(distances) + min(distances(distances > 0))).^2;
  % Method by which to find synapse:
  default_options.synapse_detection_method = 'curvature local minima';

  if ~exist('options', 'var')
    options = default_options;
  else
    options = process_options_structure(default_options, options);
  end

  
  % Things required for multiple synapse_detection_method values:

  component = single_cell_mesh;
  number_faces = size(component.faces, 1);
  number_vertices = size(component.vertices, 1);
  
  % distance_type = 'combinatorial';
  distance_type = 'distance';
  % distance_type = 'conformal';

  mesh_laplacian = compute_mesh_laplacian(component.vertices, component.faces, distance_type);
  % mesh_laplacian = -mesh_laplacian;

  
  % Compute curvature:
  [minimum_curvature_direction, maximum_curvature_direction, minimum_curvature, maximum_curvature, mean_curvature, gaussian_curvature, normal] = compute_curvature(component.vertices, component.faces, struct('curvature_smoothing', options.curvature_smoothing));
  minimum_curvature_direction = minimum_curvature_direction';
  maximum_curvature_direction = maximum_curvature_direction';
    
  % Smooth/diffuse along the mesh:
  naver = ceil(mean(std(component.vertices, 1)));
  % naver = ceil(mean(std(component.vertices, 1)) * 2);
  minimum_curvature_smoothed = perform_mesh_smoothing(component.faces, component.vertices, minimum_curvature, struct('naver', naver, 'type', distance_type));
  maximum_curvature_smoothed = perform_mesh_smoothing(component.faces, component.vertices, maximum_curvature, struct('naver', naver, 'type', distance_type));
  normal_smoothed = perform_mesh_smoothing(component.faces, component.vertices, normal, struct('naver', naver, 'type', distance_type));
  % Normalize normals:
  normal_smoothed = normal_smoothed ./ repmat(sqrt(sum(normal_smoothed.^2, 2)), 1, 3);
    

  % should_flip_normals = false;
  should_flip_normals = true;
  if should_flip_normals
    % Check if moving a small distance places this point inside the mesh (does not work for very sharp angles, but our meshes should not have those at the synapse point if they are good segmentations):
    intervertex_distances = pdist(component.vertices);
    test_vertex = 1;
    % test_point = normal_smoothed(current_synapse(1), :) .* min(intervertex_distances) * .5 + component.vertices(current_synapse(1), :);
    test_point = normal_smoothed(test_vertex, :) .* min(intervertex_distances) * .5 + component.vertices(test_vertex, :);
    if inpolyhedron(component, test_point, 'flipnormals', true)
      fprintf('>>>> Reversing normals!\n')
      normal = -normal;
      normal_smoothed = -normal_smoothed;
      % synapse_normal = -synapse_normal;
    end
  else
    fprintf('>>>> HACK, should_flip_normals = false!\n')
  end
        

  switch options.synapse_detection_method
  
  
    case 'curvature local minima'
  
  
      % Find source and target (source might be better if it were at the vertical center of mass at this XY position):
      % synapse_middle = [synapse_point(1:2), size(image_to_segment, 3) / 2]
      synapse_middle = synapse_point;
      distances = pdist2(synapse_middle, component.vertices);
      distances_2d = reshape(pdist2(synapse_middle(:, 1:2), component.vertices(:, 1:2)), [], 1);
      [~, source] = min(distances);
      [~, target] = max(distances);
      
      % Run fast marching (like the distance transform of the surface with respect to the manual synapse point and given curvature-related difficulty of travel):
      curvature_measure = minimum_curvature;
      [distances, ~, closest_points] = perform_fast_marching_mesh(component.vertices, component.faces, source, struct('W', contrast_stretch(curvature_measure)));
      coordinates = component.vertices;
      coordinates = coordinates - repmat(synapse_middle, size(coordinates, 1), 1);

      % % Robust fitting of a plane (unused now):
      % % predicted = 2; predictors = [1, 3];
      % predicted = 1; predictors = [2, 3];
      % if ischar(options.weighting_function)
        % plane_parameters = robustfit(coordinates(:, predictors), coordinates(:, predicted), options.weighting_function);
      % else
        % plane_parameters = robustfit(coordinates(:, predictors), coordinates(:, predicted), @(residuals)options.weighting_function(residuals, distances), 1.);
        % % plane_parameters = robustfit(coordinates(:, predictors), coordinates(:, predicted), @(residuals)options.weighting_function(residuals, distances_2d), 1.);
      % end

      % Try to find local optima in curvature:
      curvature_to_smooth = minimum_curvature;
      % curvature_to_smooth = gaussian_curvature;
      % Smooth/diffuse along the mesh to reduce insignificant optima:
      curvature_smoothed = perform_mesh_smoothing(component.faces, component.vertices, curvature_to_smooth, struct('naver', naver, 'type', distance_type));
      % curvature_smoothed2 = perform_mesh_smoothing(component.faces, component.vertices, curvature_to_smooth, struct('naver', ceil(mean(std(component.vertices, 1))), 'type', 'conformal'));
      % Find vertices greater than all their neighbors:
      curvature_minima = logical(process_mesh_neighborhoods(component, curvature_smoothed, @(x, y)all(x < y)));
      % curvature_minima = process_mesh_neighborhoods(component, curvature_smoothed, @(x, y)x < mean(y));
      % curvature_minima = process_mesh_neighborhoods(component, curvature_smoothed, @(x, y)x < prctile(y, 25));
      % curvature_minima = process_mesh_neighborhoods(component, curvature_smoothed, @(x, y)x / prctile(y, 25));
      
      curvature_minimum_qualities = process_mesh_neighborhoods(component, curvature_smoothed, @(x, y)prctile(y, 25) / x);
      curvature_minimum_qualities = process_mesh_neighborhoods(component, curvature_minimum_qualities, @(x, y)max(x, max(y)));
      
      % Try making quality based on alignment with the synapse point:
      mean_vertex = mean(component.vertices, 1);
      mean_direction = synapse_middle - mean_vertex;
      directions = repmat(mean_vertex, number_vertices, 1) - component.vertices;
      directions = -directions;
      directions = directions ./ repmat(sqrt(sum(directions.^2, 2)), 1, 3);
      direction_qualities = directions * mean_direction';
      direction_qualities = direction_qualities .* (direction_qualities > 0);
      direction_qualities = direction_qualities.^2;
      curvature_minimum_qualities = curvature_minimum_qualities .* direction_qualities;

      [~, synapse_vertex] = max(curvature_minimum_qualities(curvature_minima));
      curvature_minima_indices = find(curvature_minima);
      synapse_vertex = curvature_minima_indices(synapse_vertex);
      curvature_minima = curvature_smoothed == curvature_smoothed(synapse_vertex);
      
      % Segment the region within some curvature of the current vertex:
      current_synapse = synapse_vertex;
      % synapse_segmentation_curvature = curvature_smoothed;
      synapse_segmentation_curvature = maximum_curvature_smoothed;
      % synapse_segmentation_curvature = maximum_curvature;
      adjacency = logical(triangulation2adjacency(component.faces, component.vertices));
      vertex_queue = find(adjacency(synapse_vertex, :));
      % curvature_threshold = 0;
      % curvature_threshold = prctile(curvature_smoothed, 25)
      curvature_threshold = prctile(synapse_segmentation_curvature, 25);
      normal_threshold = pi / 4;
      synapse_normal = normal_smoothed(synapse_vertex, :);
      while true
        % current_synapse
        neighbors = find(max(adjacency(current_synapse, :), [], 1));
        neighbors = neighbors(~ismember(neighbors, current_synapse));
        % keyboard
        neighbor_curvatures = synapse_segmentation_curvature(neighbors);
        % pause
        % neighbors_to_add = neighbors(neighbor_curvatures < curvature_threshold);
        % % Decrease the threshold as the patch grows:
        % current_curvature_threshold = curvature_threshold * (number_vertices - length(neighbors)) / number_vertices;
        % neighbors_to_add = neighbors(neighbor_curvatures < curvature_threshold);
        % scalar_products = normal_smoothed(neighbors, :) * synapse_normal';
        scalar_products = max(min(normal_smoothed(neighbors, :) * synapse_normal', 1), -1);
        neighbor_angles = acos(scalar_products);
        neighbors_to_add = (neighbor_curvatures <= curvature_threshold) & (neighbor_angles <= normal_threshold);
        novel_neighbors = neighbors(neighbors_to_add);
        if length(novel_neighbors) == 0
          break
        end
        % Just use the best new neighbor:
        neighbor_curvatures = neighbor_curvatures(neighbors_to_add);
        scalar_products = scalar_products(neighbors_to_add);
        neighbor_angles = neighbor_angles(neighbors_to_add);
        [~, neighbor_curvatures_order] = sort(neighbor_curvatures);
        [~, neighbor_angles_order] = sort(neighbor_angles);
        % Minimum angle to synapse:
        % [~, novel_neighbor_index] = min(neighbor_angles);
        % novel_neighbors = novel_neighbors(novel_neighbor_index);
        % Heuristic minimum angle to synapse and minimum maximum curvature:
        [~, novel_neighbor_index] = min(neighbor_curvatures_order + neighbor_angles_order);
        novel_neighbors = novel_neighbors(novel_neighbor_index);
        
        % Add the new neighbors and update the normal:
        current_synapse = [current_synapse, novel_neighbors];
        synapse_normal = mean(normal_smoothed(current_synapse, :), 1);
        synapse_normal = synapse_normal ./ norm(synapse_normal);
      end
      in_synapse = false(number_vertices, 1);
      in_synapse(current_synapse) = true;
      synapse_center = mean(component.vertices(in_synapse, :), 1);
      

      % Fit a fully 3D plane:
      synapse_centered_vertices = (component.vertices - repmat(synapse_center, number_vertices, 1));
      synapse_centered_distances = synapse_centered_vertices * synapse_normal';
      % flattened_vertices = synapse_centered_vertices - repmat(synapse_centered_distances, 1, 3);
      flattened_vertices = synapse_centered_vertices - synapse_centered_distances * synapse_normal;
      flattened_bases = flattened_vertices(1, :) - flattened_vertices(2, :);
      flattened_bases = flattened_bases ./ norm(flattened_bases);
      flattened_bases = [flattened_bases; cross(flattened_bases, synapse_normal)];
      % whos flattened_vertices flattened_bases
      flattened_vertices = flattened_vertices * flattened_bases';
      
      % New radius estimate (maximum distance in the synapse plane across vertices with positive minimum curvature and within one third of the maximum distance of a vertex from the synapse plane):
      high_curvature_vertices = (minimum_curvature_smoothed > 0) & (synapse_centered_distances < max(synapse_centered_distances) / 3);
      synapse_radius = max(sqrt(sum(flattened_vertices(high_curvature_vertices, :).^2, 2)));
      
      % keyboard
  

    case 'synapse direction'
    
    
      % % Use the direction from the mean vertex to the synapse as the orientation:
      % mean_vertex = mean(single_cell_mesh.vertices, 1);
      % synapse_normal = synapse_point - mean_vertex;
      
      % Use the direction from the centroid to the synapse as the orientation:
      mesh_integer_ranges = [floor(min(single_cell_mesh.vertices, [], 1)); ceil(max(single_cell_mesh.vertices, [], 1))];
      current_x_range = mesh_integer_ranges(1, 1):mesh_integer_ranges(2, 1);
      current_y_range = mesh_integer_ranges(1, 2):mesh_integer_ranges(2, 2);
      current_z_range = mesh_integer_ranges(1, 3):mesh_integer_ranges(2, 3);
      segmentation_image = voxelise(current_x_range, current_y_range, current_z_range, single_cell_mesh);
      % segmentation_image = permute(segmentation_image, [2, 1, 3]) > 0;
      segmentation_image = permute(segmentation_image, [2, 1, 3]);
      % volume = sum(segmentation_image(:))
      [x, y, z] = meshgrid(current_x_range, current_y_range, current_z_range);
      % centroid = [sum(x(segmentation_image)) ./ volume, sum(y(segmentation_image)) ./ volume, sum(z(segmentation_image)) ./ volume]
      centroid = mean([x(segmentation_image), y(segmentation_image), z(segmentation_image)], 1);
      synapse_normal = synapse_point - centroid;

      synapse_normal = synapse_normal ./ norm(synapse_normal);
      
      distances_2d = reshape(pdist2(synapse_point(:, 1:2), single_cell_mesh.vertices(:, 1:2)), [], 1);
      % nearby_vertices_2d = distances_2d <= prctile(distances_2d, 50);
      
      % Project points onto this plane:
      % We don't know synapse_center yet:
      % synapse_centered_vertices = (component.vertices - repmat(synapse_center, number_vertices, 1));
      synapse_centered_vertices = component.vertices;
      synapse_centered_distances = synapse_centered_vertices * synapse_normal';
      flattened_vertices = synapse_centered_vertices - synapse_centered_distances * synapse_normal;
      flattened_bases = flattened_vertices(1, :) - flattened_vertices(2, :);
      flattened_bases = flattened_bases ./ norm(flattened_bases);
      flattened_bases = [flattened_bases; cross(flattened_bases, synapse_normal)];
      flattened_vertices_2d = flattened_vertices * flattened_bases';
      
      % nearby_vertices_percentile = 50;
      % This is no longer a hack, it is an OK idea for most cells:
      % fprintf('>>>> HACK, nearby_vertices_percentile = 0!\n'), nearby_vertices_percentile = 0;      
      nearby_vertices_percentile = 0;      
      % nearby_vertices = distances_2d <= prctile(distances_2d, nearby_vertices_percentile);
      nearby_vertices = synapse_centered_distances >= prctile(synapse_centered_distances, nearby_vertices_percentile);
      
      % fprintf('>>>> HACK, X and Y switched!\n'), flattened_vertices_2d = flattened_vertices_2d(:, [2, 1]);
      
      % Compute the bounding ellipsoid of the nearest points:
      [lowner_ellipse, lowner_ellipse_center] = lowner(flattened_vertices_2d(nearby_vertices, :)', 1e-3);
      [V D]=eig(lowner_ellipse);
      t=linspace(0,2*pi);
      e=repmat(lowner_ellipse_center,size(t))+V*diag(1./sqrt(diag(D)))*[cos(t);sin(t)];

      % The synapse center is the ellipse center at the 95th percentile distance along the normal of nearby vertices:
      synapse_center = (flattened_bases' * lowner_ellipse_center)' + synapse_normal * prctile(synapse_centered_distances(nearby_vertices), 95);
      
      [~, major_axis_index] = min(diag(D));
      lowner_ellipse_axes = zeros(2, 2); 
      lowner_ellipse_axes(major_axis_index, :) = (V * diag(1./sqrt(diag(D))) * [1; 0])';
      lowner_ellipse_axes(3 - major_axis_index, :) = (V * diag(1./sqrt(diag(D))) * [0; 1])';
      
      synapse_radius = (norm(lowner_ellipse_axes(major_axis_index, :)) + norm(lowner_ellipse_axes(3 - major_axis_index, :))) / 2;

      % synapse_center, synapse_normal, synapse_radius

      % % Show synapse center and normal:
      % % scatter3(single_cell_mesh.vertices(:, 1), single_cell_mesh.vertices(:, 2), single_cell_mesh.vertices(:, 3), 'r.'), hold on, scatter3(synapse_center(1), synapse_center(2), synapse_center(3), 'gx'), quiver3(synapse_center(1), synapse_center(2), synapse_center(3), synapse_normal(1), synapse_normal(2), synapse_normal(3), 5, 'b'), hold off, axis equal, pause
      % clf
      % patch(single_cell_mesh, 'EdgeColor', 'none', 'FaceColor', [1, .4, .1])
      % hold on
      % scatter3(synapse_center(1), synapse_center(2), synapse_center(3), 'gx')
      % quiver3(synapse_center(1), synapse_center(2), synapse_center(3), synapse_normal(1), synapse_normal(2), synapse_normal(3), 5, 'b')
      % hold off
      % axis equal
      % light_handle = lightangle(0, 60); set(light_handle, 'Color', .75 * ones(1, 3));
      % light_handle = lightangle(150, -30); set(light_handle, 'Color', .5 * ones(1, 3));
      % lighting phong
      % pause
      
      % % Show points parallel to the synapse plane along with arrow pointing along minor axis:
      % plot(flattened_vertices_2d(nearby_vertices, 1), flattened_vertices_2d(nearby_vertices, 2),'+',e(1,:),e(2,:),'-')
      % hold on
      % % major_minor_points = repmat(lowner_ellipse_center', 4, 1) + [synapse_radius(:, 1)'; synapse_radius2(:, 1)'; synapse_radius(:, 2)'; synapse_radius2(:, 2)']
      % % major_minor_points = repmat(lowner_ellipse_center', 3, 1) + [synapse_radius(:, 1)'; synapse_radius2(:, 1)'; synapse_radius(:, 2)']
      % major_minor_points = repmat(lowner_ellipse_center', 3, 1) + [lowner_ellipse_axes(1, :); lowner_ellipse_axes(2, :); -lowner_ellipse_axes(1, :)]
      % plot(major_minor_points(:, 1), major_minor_points(:, 2), 'r-')
      % hold off
      % axis equal
      % pause
      
      % keyboard
  
    otherwise
      error('Unknown synapse_detection_method value ''%s''', options.synapse_detection_method)
  end
  
  
  % plane_parameters = [synapse_center, synapse_normal, synapse_radius];
  plane_parameters = struct();
  plane_parameters.synapse_center = synapse_center;
  plane_parameters.synapse_normal = synapse_normal;
  % Should eventually use an ellipse instead of a circle:
  % use_computed_synapse_radius = false;
  % % use_computed_synapse_radius = true;
  % if use_computed_synapse_radius
    plane_parameters.synapse_radius = synapse_radius;
  % else
    % plane_parameters.synapse_radius = 1;
    % fprintf('>>>> HACK, use_computed_synapse_radius = false!\n')
  % end
  % plane_parameters
  
  if isnan(synapse_radius)
    % keyboard
    error('synapse_radius is nan!')
  end
  

  % Plot the segmentation surface and the split:
  if options.debug_plot
    clf
    hold on
    light_handle = lightangle(0, 60); set(light_handle, 'Color', .75 * ones(1, 3));
    light_handle = lightangle(150, -30); set(light_handle, 'Color', .5 * ones(1, 3));

    component2 = component; component2.vertices(:, 1) = component2.vertices(:, 1) + range(component2.vertices(:, 1)) * 1.5;
    component3 = component; component3.vertices(:, 1) = component3.vertices(:, 1) + range(component3.vertices(:, 1)) * 1.5 * 2;
    if exist('curvature_smoothed', 'var')
      first_mesh_plot_handle = patch(component, 'CData', contrast_stretch(curvature_smoothed), 'EdgeColor', 'none', 'FaceColor', 'interp');
      patch(component2, 'CData', contrast_stretch(curvature_minima), 'EdgeColor', 'none', 'FaceColor', 'interp')
      patch(component3, 'CData', contrast_stretch(in_synapse), 'EdgeColor', 'none', 'FaceColor', 'interp')
    else
      first_mesh_plot_handle = patch(component, 'EdgeColor', 'none', 'FaceColor', [1, .4, .1]);
    end
    
    z_spacing = (max(component.vertices(:, 3)) - min(component.vertices(:, 3))) * 1.5;

    lighting phong

    plot3(synapse_point(1) * [1, 1], synapse_point(2) * [1, 1], synapse_point(3) + z_spacing * [-.5, .5], 'm-')
    
    if exist('curvature_smoothed', 'var')
      quiver_positions = repmat(synapse_center, 3, 1) + [(0:2)' * range(component2.vertices(:, 1)) * 1.5, zeros(3, 2)];
      quiver3(quiver_positions(:, 1), quiver_positions(:, 2), quiver_positions(:, 3), synapse_normal(1) * ones(3, 1), synapse_normal(2) * ones(3, 1), synapse_normal(3) * ones(3, 1), 'g')
    else
      quiver_positions = synapse_center;
      % quiver3(quiver_positions(:, 1), quiver_positions(:, 2), quiver_positions(:, 3), synapse_normal(1), synapse_normal(2), synapse_normal(3), 'm')
      % quiver3(quiver_positions(:, 1), quiver_positions(:, 2), quiver_positions(:, 3), synapse_normal(1), synapse_normal(2), synapse_normal(3), 10, 'm')
      synapse_normal_stretched = synapse_normal * 20;
      % synapse_normal_stretched = synapse_normal * 100;
      quiver3(quiver_positions(:, 1), quiver_positions(:, 2), quiver_positions(:, 3), synapse_normal_stretched(1), synapse_normal_stretched(2), synapse_normal_stretched(3), 0, 'g')
    end
    
    % if exist('centroid', 'var')
      % hold on
      % scatter3(centroid(1), centroid(2), centroid(3), 'bo')
      % hold off
      % % set(first_mesh_plot_handle, 'EdgeColor', [.1, .4, 1], 'FaceColor', 'none', 'EdgeLighting', 'phong')
      % set(first_mesh_plot_handle, 'EdgeColor', [.3, 1, .2], 'FaceColor', 'none', 'EdgeLighting', 'phong')
    % end

    % view(3)
    [az, el, r] = cart2sph(synapse_normal(1), synapse_normal(2), synapse_normal(3));
    % [az, el, r] = cart2sph(synapse_normal(1), synapse_normal(3), synapse_normal(2))
    % view(az * 180 / pi + 90, el * 180 / pi)
    % This way we can actually see the synapse normal:
    view(az * 180 / pi + 90 + 30, el * 180 / pi)
    axis equal
    % zlim(synapse_point(3) + z_spacing * [-.5, .5])

    hold off
  end
  
  

  
% % Visualize:
% faces_to_remove = vertices_to_remove(synapse_mesh.faces(:, 1)) | vertices_to_remove(synapse_mesh.faces(:, 2)) | vertices_to_remove(synapse_mesh.faces(:, 3));
% number_faces_to_remove = sum(faces_to_remove);
% number_vertices_to_remove = sum(vertices_to_remove);
% synapse_mesh.faces = synapse_mesh.faces(~faces_to_remove, :);
% clf
% patch(synapse_mesh, 'FaceColor', [1, .4, .1], 'EdgeColor', 'interp', 'CData', synapse_curvature_transformed, 'SpecularStrength', 0), camlight, lighting phong, axis equal, view(0, 0)

