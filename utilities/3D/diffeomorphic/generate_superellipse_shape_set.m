function [shapes, shape_parameters] = generate_superellipse_shape_set(options)
% Create a set of binary (or optionally antialiased) images of superelliptical shapes for testing shape space-related computations.


  default_options.image_width = 64;
  default_options.number_images = 64;
  % Generate the image at this many times the width, then downsample to antialias:
  default_options.antialiasing_factor = 3;
  default_options.minimum_relative_semidiameter = 1 / 4;
  default_options.maximum_relative_semidiameter = 2 / 3;
  default_options.minimum_exponent = 3 / 4;
  default_options.maximum_exponent = 8;
  default_options.interpolation_factor_powers = [1, 2];
  default_options.generate_cycle = false;
  default_options.generate_randomly = false;
  default_options.debug = false;
  
  if ~exist('options', 'var')
    options = default_options; 
  else
    options = process_options_structure(default_options, options); 
  end
  
  

  % Generate some shapes (just bars for now):
  shapes = cell(options.number_images, 1);
  shape_parameters = cell(options.number_images, 1);

  % Create ellipses:
  [x, y] = meshgrid(1:(options.image_width * options.antialiasing_factor), 1:(options.image_width * options.antialiasing_factor));
  a = mean([options.minimum_relative_semidiameter, options.maximum_relative_semidiameter]) * options.image_width;
  for shape_index = 1:options.number_images
    if options.generate_cycle
      % One-parameter cyclic family of superellipses:
      if options.generate_randomly
        interpolation_factors = rand() * 2 * pi;
      else
        interpolation_factors = (shape_index - 1) / options.number_images * 2 * pi;
      end
      interpolation_factors = [cos(interpolation_factors), sin(interpolation_factors)] * .5 + .5;
    else
      % Two-parameter family of superellipses:
      if options.generate_randomly
        interpolation_factors = rand(1, 2);
      else
        interpolation_factors = [mod(shape_index - 1, sqrt(options.number_images)) / (sqrt(options.number_images) - 1), floor((shape_index - 1) / sqrt(options.number_images)) / (sqrt(options.number_images) - 1)];
      end
    end
    
    interpolation_factors = interpolation_factors.^options.interpolation_factor_powers;
    % interpolation_factors
      
    b = (options.minimum_relative_semidiameter .* (1 - interpolation_factors(1)) + options.maximum_relative_semidiameter .* interpolation_factors(1)) * options.image_width;
    c = options.minimum_exponent .* (1 - interpolation_factors(2)) + options.maximum_exponent .* interpolation_factors(2);
    
    % [a, b, c]
    
    shapes{shape_index} = abs((x - (options.image_width * options.antialiasing_factor)/2) / a).^c + abs((y - (options.image_width * options.antialiasing_factor)/2) / b).^c - 1 <= 0;
    shapes{shape_index} = imresize(double(shapes{shape_index}), [1, 1] * options.image_width, 'bilinear');
    shape_parameters{shape_index} = [b, c];
    
    if options.debug
      imshow(shapes{shape_index}), pause
    end
  end

  if options.debug
    error
  end
  
  shape_parameters = cell2mat(shape_parameters);


  

end
