function [value] = image_central_moment(image, p, q, center)
  %fprintf('image_central_moment\n')
  %center
  [h, w] = size(image);
  if ~exist('center', 'var') || isempty(center)
    m00 = image_moment(image, 0, 0);
    m10 = image_moment(image, 1, 0);
    m01 = image_moment(image, 0, 1);
    center = [m10 / m00, m01 / m00];
  else
    center = center ./ [w, h];
  end
  %center
  
  value = sum(((1:h) ./ h - center(2)).^q * image, 1);
  value = value * (((1:w) ./ w - center(1)).^p)';
