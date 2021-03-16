function [centroid] = image_centroid(given_image)
  % Centroid returned in XYZ order.
  %
  % Tests:
  % image_centroid(ones(5, 10, 15))
  % image_centroid(rand(5, 10, 15))
  
  number_dimensions = ndims(given_image);
  centroid = zeros(1, number_dimensions);
  total_intensity = sum(given_image(:));
  given_image = given_image ./ total_intensity;
  for dimension_index = 1:number_dimensions
    dimension_range = (1:size(given_image, dimension_index))';
    dimension_range = permute(dimension_range, circshift((1:number_dimensions)', dimension_index - 1));
    repeat_size = size(given_image);
    repeat_size(dimension_index) = 1;
    a = repmat(dimension_range, repeat_size);
    centroid(dimension_index) = sum(given_image(:) .* a(:));
  end
  
  centroid = centroid([2, 1, 3:end]);
  
end
