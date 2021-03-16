function [r_image] = contrast_stretch(given_image)
  r_image = (given_image - min(given_image(:))) ./ (max(given_image(:)) - min(given_image(:)));
  %r_image = (double(given_image) - min(given_image(:))) ./ (max(given_image(:)) - min(given_image(:)));
