function [b] = reshape_contrast(a, r, is_vector_field)
  % Convenience function to apply contrast_stretch to a before applying reshape_2d. 
  if ~exist('r', 'var') || isempty(r)
    r = 1;
  end
  if ~exist('is_vector_field', 'var') || isempty(is_vector_field)
    is_vector_field = false;
  end
  if is_vector_field
    maximum_absolute_value = max(abs(a(:)));
    for index = 1:3
      b(:, :, index) = reshape_2d(a(:, :, :, index), r);
      b(:, :, index) = contrast_stretch_about_zero(b(:, :, index));
    end
  else
    b = a;
    b = contrast_stretch(b);
    %D. Sullivan 12/15/24 - when synthesizing 2D images, no need to do
    %this.
    if ndims(b)==3
        b = reshape_2d(b, r);
    end
  end
