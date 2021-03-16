function [result_image] = image_resize_nd(given_image, scale, method, use_scale_not_size)
  image_dimensionality = max(ndims(given_image), numel(scale)); 
  if isscalar(scale)
    scale = scale * ones(1, image_dimensionality); 
  end
  if ~exist('method', 'var')
    method = 'bicubic';
  end
  if ~exist('use_scale_not_size', 'var')
    use_scale_not_size = false;
  end
  
  [discard, dimension_order] = sort(scale);
  if use_scale_not_size
    for dimension_index = dimension_order
      % Reorder the dimensions so that scaling along the first
      % dimension/vertical axis of the given_image is resizing along that
      % dimension after returning the given_image to the original order:
      current_dimension_order = [dimension_index, 1:(dimension_index - 1), (dimension_index + 1):image_dimensionality]; 
      given_image = permute(given_image, current_dimension_order); 
      image_size = size(given_image);
      novel_image_sizes = num2cell(image_size);
      novel_image_sizes{1} = [];
      given_image = reshape(given_image, image_size(1), image_size(2), []);
      % given_image = imresize(given_image, scale(dimension_index), method); 
      % given_image = imresize(given_image, [size(given_image, 1), image_size(2)], method); 
      % given_image = imresize(given_image, scale(dimension_index), method, 'Scale', [scale(dimension_index), 1]); 
      given_image = imresize(given_image, method, 'Scale', [scale(dimension_index), 1]); 
      given_image = reshape(given_image, novel_image_sizes{:});
      given_image = ipermute(given_image, current_dimension_order); 
      if false
        % Debug info:
        fprintf('dimension_index %d, scale %f\n', dimension_index, scale(dimension_index))
      end
    end
  else
    for dimension_index = dimension_order
      % Reorder the dimensions so that scaling along the first
      % dimension/vertical axis of the given_image is resizing along that
      % dimension after returning the given_image to the original order:
      current_dimension_order = [dimension_index, 1:(dimension_index - 1), (dimension_index + 1):image_dimensionality]; 
      given_image = permute(given_image, current_dimension_order); 
      image_size = size(given_image);
      novel_image_size = image_size;
      novel_image_size(1:2) = ceil(novel_image_size(1:2) .* [scale(dimension_index), 1]);
      given_image = reshape(given_image, image_size(1), image_size(2), []);
      given_image = imresize(given_image, novel_image_size(1:2), method); 
      given_image = reshape(given_image, novel_image_size);
      given_image = ipermute(given_image, current_dimension_order); 
    end
  end
  
  result_image = given_image;
  
end
  