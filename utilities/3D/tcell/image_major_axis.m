function [value] = image_major_axis(image, center)
  if ~exist('center', 'var')
    center = [];
  end
  
  mu11 = image_central_moment(image, 1, 1, center);
  mu20 = image_central_moment(image, 2, 0, center);
  mu02 = image_central_moment(image, 0, 2, center);
% $$$   value = .5 * atan2(2 * mu11 / mu00, mu20 / mu00 - mu02 / mu00);
  value = .5 * atan2(2 * mu11, mu20 - mu02);
    
  
  