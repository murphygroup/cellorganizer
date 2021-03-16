function [] = illustrate_shape_space(options)
% Plot training shapes of a shape space model.
% 
% Options you must specify for good results:
% image_function
% voxel_size
% 
% Options you should specify for completeness:
% 
% 
% Other options can be seen in the default_options structure at the beginning of this function.
% 
% Dependencies from File Exchange: isocontour_version2, pmkmp
% 
% History:
% 2013-02-14 tebuck: Created.
% 2013-06-28 tebuck: Copying over more interesting code from compute_hela3d_willmore_2d_shape_space.m.
% 
% To do:
% Allow using meshes instead of images.


  % We require certain options:
  required_options = {'image_function', 'voxel_size', 'shape_space', 'number_images'};
  for index = 1:length(required_options)
    required_option = required_options{index};
    if ~isfield(options, required_option)
      error('Option %s is required', required_option)
    end
  end


  first_image = options.image_function(1);

  [M, N, P] = size(first_image);
  [X_0, Y_0, Z_0] = meshgrid(1:N, 1:M, 1:P);
  identity_map = {X_0, Y_0, Z_0}; 
  clear X_0 Y_0 Z_0
  
  window_size = [M, N, P];

  
  default_options.marker_size_scale = 2;
  % Should add an option for showing isosurfaces instead of isocontours:
  default_options.shape_projection_function = @(x)max(x, [], 3);
  % Levels in the image that correspond to different shapes, with the lowest assumed background, e.g., [0, 1] for an image where closer to 0 is background and closer to 1 is foreground:
  default_options.shape_levels = [0, 1, 2];
  default_options.individual_marker_size_influence = 1;
  default_options.use_outlines = false;
  default_options.use_third_coordinate_for_color = false;
  default_options.shape_space_inlier_quantile = 1;
  default_options.shape_parameters = {};
  

  options = process_options_structure(default_options, options); 
  
  
  
  positions = options.shape_space.positions;
  number_images = options.number_images;
%   distances = options.shape_space.distances;
%   
  
  distances_2d = pdist(positions(:, 1:2));
  

  inliers = mahal(positions, positions);
  inliers = inliers <= quantile(inliers, options.shape_space_inlier_quantile);
  positions_without_outliers = positions;
  positions_without_outliers(~inliers, :) = nan;
  % inliers
  
  space_limits_x = [min(positions_without_outliers(inliers, 1)), max(positions_without_outliers(inliers, 1))];
  space_limits_y = [min(positions_without_outliers(inliers, 2)), max(positions_without_outliers(inliers, 2))];
    





  % % marker_size = mean(std(positions)) * options.marker_size_scale * sqrt(50 / number_images); 
  % % marker_size = mean(std(positions(:, 1:2))) * options.marker_size_scale * prctile(distances(distances > 0), 50);
  % % marker_size = options.marker_size_scale * prctile(distances(distances > 0), 50);
  % marker_size = options.marker_size_scale * prctile(distances_2d(distances_2d > 0), 50);
  % marker_size
  % % cellular_colors = hsv2rgb([rand(number_images, 1), rand(number_images, 1) * .6 + .4, .5 * ones(number_images, 1)]);
  % cellular_colors = hsv2rgb([rand(number_images, 1), rand(number_images, 1) * .6 + .4, 1 * ones(number_images, 1)]);

  
  marker_size = mean(iqr(positions(:, 1:2))) * options.marker_size_scale * sqrt(50 / number_images); 
  
  if options.individual_marker_size_influence > 0
    % Scale shapes based on nearest neighbor distances:
    % individual_marker_sizes = sort(distances);
    % 2D:
    individual_marker_sizes = sort(squareform(distances_2d), 2);
    % % 3D:
    % individual_marker_sizes = sort(squareform(pdist(positions(:, 1:3))), 2);
    individual_marker_sizes = individual_marker_sizes(:, 2);
    individual_marker_sizes = options.individual_marker_size_influence * individual_marker_sizes + (1 - options.individual_marker_size_influence) * marker_size;
  
    % Could use a hyperbolic projection or log of radial distance from shape space center... Or do KDE and solve for distorted space that best expands interesting areas?
  else
    individual_marker_sizes = ones(number_images, 1) * options.marker_size_scale * prctile(distances_2d(distances_2d > 0), 50);
  end
  
  % individual_marker_sizes

  % Choose individual colors:
  if (size(positions, 2) > 2 && options.use_third_coordinate_for_color) || ~isempty(options.shape_parameters)
    % Color indicates third coordinate:
    % uniform_colors = pmkmp(1024);
    uniform_colors = pmkmp(1024, 'IsoL');
    if (size(positions, 2) > 2 && options.use_third_coordinate_for_color)
      cellular_colors = interp1((0:1023) ./ 1023 .* range(positions(inliers, 3)) + min(positions(inliers, 3)), uniform_colors, positions(:, 3), 'linear', 'extrap');
    else
      disp('Using parameters for colors!')
      cellular_colors = interp1((0:1023) ./ 1023 .* range(options.shape_parameters(inliers, 1)) + min(options.shape_parameters(inliers, 1)), uniform_colors, options.shape_parameters(inliers, 1), 'linear', 'extrap');
      second_interpolation_factors = (options.shape_parameters(inliers, 2) - min(options.shape_parameters(inliers, 2))) ./ range(options.shape_parameters(inliers, 2));
      % second_interpolation_factors
      second_interpolation_factors_repeated = repmat(second_interpolation_factors, [1, 3]);
      % Increasing shape_parameters(:, 2) saturates the color
      desaturation_effect = .5;
      % desaturation_effect = 1;
      % cellular_colors = cellular_colors .* (second_interpolation_factors_repeated * .5 + .5) + (1 - second_interpolation_factors_repeated) .* .5;
      cellular_colors = cellular_colors .* (second_interpolation_factors_repeated * (desaturation_effect) + (1 - desaturation_effect)) + (1 - second_interpolation_factors_repeated) .* desaturation_effect;
    end
  else
    cellular_colors = hsv2rgb([rand(number_images, 1), rand(number_images, 1) * .6 + .4, 1 * ones(number_images, 1)]);
  end

  number_shape_levels = length(options.shape_levels);
  
  % figure
  % set(gcf, 'Color', 'k');
  % set(gca, 'Color', 'none');
  % % set(gca, 'Color', 'k');
  % set(gca, 'XColor', 'w');
  % set(gca, 'YColor', 'w');
  % set(gca, 'ZColor', 'w');
  % hold on
  
  % for n = 1:options.number_images
    % datum_shape = options.shape_projection_function(options.image_function(n));
    % for shape_index = 1:number_shape_levels - 1
      % current_threshold = mean(options.shape_levels(shape_index:shape_index + 1));
      % [lines, vertices, objects] = isocontour(datum_shape, current_threshold);
      % datum_polygon = vertices(objects{1}, :);
      % datum_polygon = datum_polygon(:, [2, 1]);
      % if options.use_outlines
        % patch((datum_polygon(:, 1) / N - .5) * marker_size + positions(n, 1), (datum_polygon(:, 2) / M - .5) * marker_size + positions(n, 2), cellular_colors(n, :) * (shape_index / (number_shape_levels - 1)), 'FaceColor', 'none', 'EdgeColor', cellular_colors(n, :) * (shape_index / (number_shape_levels - 1)));
      % else
        % fill((datum_polygon(:, 1) / N - .5) * marker_size + positions(n, 1), (datum_polygon(:, 2) / M - .5) * marker_size + positions(n, 2), cellular_colors(n, :) * (shape_index / (number_shape_levels - 1)), 'EdgeColor', 'none');
      % end
    % end
  % end
  
  % % set(gca, 'YDir', 'reverse')
  % axis image
  % hold off
  % % colormap(gray)


  clf
  % set(gcf, 'Visible', 'off');
  % screen_dpi = get(0, 'ScreenPixelsPerInch');
  % set(gcf, 'Position', [1, 1, 0, 0] + [0, 0, 6, 4] * screen_dpi);
  % set(gcf, 'PaperPositionMode', 'auto');
  % set(gcf, 'Color', 'k');
  % whitebg
  % set(gca, 'Color', [0, 0, 0]);
  % set(gca, 'Color', 'none');
  % set(gca, 'XColor', 'w');
  % set(gca, 'YColor', 'w');
  % set(gca, 'ZColor', 'w');
  % set(subplot_handle, 'Color', 'none');
  % set(subplot_handle, 'XColor', 'w');
  % set(subplot_handle, 'YColor', 'w');
  % set(subplot_handle, 'ZColor', 'w');
  hold on
  
  for image_index = 1:number_images
    % if ~all(isfinite(positions_without_outliers(image_index, :)))
    if ~inliers(image_index, :)
      continue
    end

    datum_shape = options.shape_projection_function(options.image_function(image_index));
    for shape_index = 1:number_shape_levels - 1
      current_threshold = mean(options.shape_levels(shape_index:shape_index + 1));
      [lines, vertices, objects] = isocontour(datum_shape, current_threshold);
      datum_polygon = vertices(objects{1}, :);
      datum_polygon = datum_polygon(:, [2, 1]);
      
      % Draw the shape:
      if options.use_outlines
        fill_extra_arguments = {'EdgeColor', cellular_colors(image_index, :), 'FaceColor', 'none'};
      else
        fill_extra_arguments = {};
      end
      fill((datum_polygon(:, 1) / N - .5) * individual_marker_sizes(image_index) + positions(image_index, 1), (datum_polygon(:, 2) / M - .5) * individual_marker_sizes(image_index) + positions(image_index, 2), cellular_colors(image_index, :), fill_extra_arguments{:})
      % patch((datum_polygon(:, 1) / N - .5) * individual_marker_sizes(image_index) + positions(image_index, 1), (datum_polygon(:, 2) / M - .5) * individual_marker_sizes(image_index) + positions(image_index, 2), cellular_colors(image_index, :))
      
      % plot(positions(image_index, 1), positions(image_index, 2), 'x')
    end

  end
  
  % set(gca, 'YDir', 'reverse')
  axis image
  hold off
  % if options.use_third_coordinate_for_color
  if (size(positions, 2) > 2 && options.use_third_coordinate_for_color) || ~isempty(options.shape_parameters)
    % colormap(gray)
    colormap(pmkmp)
    colorbar_axes = colorbar;
    if (size(positions, 2) > 2 && options.use_third_coordinate_for_color)
      colormap_label = 'Component 3';
    else
      colormap_label = {'Hue: first shape parameter', 'Saturation: second shape parameter'};
    end
    ylabel(colorbar_axes, colormap_label, 'Interpreter', 'none')
    % ylabel(colorbar_axes, 'Curvature')
    % caxis([min(positions(inliers, 3)), max(positions(inliers, 3))])
  end

  % xlim(space_limits_x);
  % ylim(space_limits_y);
  
  xlabel('Component 1')
  ylabel('Component 2')
  
  if diff(space_limits_y) > diff(space_limits_x)
    view(-90, 90)
    set(gca, 'YDir', 'reverse')
  end


  

end
