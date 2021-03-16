function [value] = image_moment(image, p, q)
  [h, w] = size(image);
  value = sum(((1:h) ./ h).^q * image, 1);
  value = value * (((1:w) ./ w).^p)';
