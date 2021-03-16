function [r_image] = contrast_stretch(image)
  % Linearly transform values such that the minimum value becomes zero and the maximum value one.
  r_image = (image - min(image(:))) ./ (max(image(:)) - min(image(:)));
  %r_image = (double(image) - min(image(:))) ./ (max(image(:)) - min(image(:)));
