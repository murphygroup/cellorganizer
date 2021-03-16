function [fx, fy, fz] = ...
      smooth_gradient(...
        image, filter_type, use_miter, use_padding)
  % Compute gradients using a variety of derivative filters.
  %
  % Requires "Steerable Filters" from Mathworks Exchange in some
  % cases.
  if (~exist('filter_type', 'var'))
    %filter_type = 'kroon';
    filter_type = 'gradient';
  end
  if (~exist('use_miter', 'var'))
    use_miter = false;
  end
  if (~exist('use_padding', 'var'))
    use_padding = false;
  end
  
  gradient_computed = false;

  miter_radius = 2;
  
  % Backwards because it's convolution, not correlation:
  gradient_filter = [1, 0, -1];
  gradient_filter = gradient_filter * .5;
  smoothing_filter = [1];
  % Use Matlab's gradient function:
  fx = []; 
  fy = []; 
  fz = []; 
  if (strcmpi(filter_type, 'gradient'))
    % if size(image, 3) == 1
      % [fx, fy] = gradient(image);
      % fz = zeros(size(fx)); 
    % else
      % [fx, fy, fz] = gradient(image);
    % end
    % gradient_computed = true;
  % Central difference, for debugging, no edge case handling:
  elseif (strcmpi(filter_type, 'central'))
    smoothing_filter = [0, 1, 0];
  % Sobel:
  elseif (strcmpi(filter_type, 'sobel'))
    smoothing_filter = [1, 2, 1];
  % Scharr 3 3D ("Optimal Filters for Extended Optical Flow", Table 3):
  elseif (strcmpi(filter_type, 'scharr'))
    %Old, from Kroon paper:
    %smoothing_filter = [3, 10, 3];
    % From Scharr's paper:
    smoothing_filter = [0.1837, 0.6326, 0.1837];
  % Scharr 5 3D:
  elseif (strcmpi(filter_type, 'scharr5'))
    gradient_filter = [0.0836, 0.3327, 0, -0.3327, -0.0836];
    % Already normalized.
    %% Normalization from same paper:
    %gradient_filter = gradient_filter .* (1. / ((-2:2) * gradient_filter'));
    smoothing_filter = [0.0233, 0.2415, 0.4704, 0.2415, 0.0233];
    miter_radius = 3;
  % Kroon (http://www.k-zone.nl/Kroon_DerivativePaper.pdf):
  elseif (strcmpi(filter_type, 'kroon'))
    smoothing_filter = [17, 61, 17];
  % Should add decomposition of Kroon's 5x5, others...
  % Sobel 5x5 (Scharr 2007, "Optimal Filters for Extended Optical Flow"):
  elseif (strcmpi(filter_type, 'sobel5'))
    gradient_filter = [-1, 8, 0, -8, 1] ./ 12.;
    smoothing_filter = [1, 4, 6, 4, 1];
    miter_radius = 3;
  % Sobel 7x7?
  % Canny-Sobel (Canny-sized Gaussian's gradient):
  elseif (strcmpi(filter_type, 'cannysobel'))
    smoothing_filter = fspecial('gaussian', [5, 1], 1.4);
    miter_radius = 3;
  % Canny (just Gaussian-smoothed):
  elseif (strcmpi(filter_type, 'canny'))
    error('Not implemented as an actual Canny filter, and should use smoothing_filter like other options')
      % image = smooth3(image, 'gaussian', 5, 1.4);
    % if size(image, 3) == 1
      % [fx, fy] = gradient(image);
      % fz = zeros(size(fx)); 
    % else
      % [fx, fy, fz] = gradient(image);
    % end
    % gradient_computed = true;
    % miter_radius = 3;
    % %return
  end
  %nvm, not 3d:
  % Steerable filters?
  
  pad_size = round(length(smoothing_filter) / 2 - .5) + round(length(gradient_filter) / 2 - .5);
  if use_padding
    image = padarray(image, ones(1, 3) * pad_size, 'replicate');
  end
  
  if ~gradient_computed
    % Make it 2D:
    smoothing_filter = smoothing_filter' * smoothing_filter; 
    smoothing_filter = smoothing_filter ./ sum(smoothing_filter(:));
    % fx:
    x_filter = convn(...
      permute(gradient_filter, circshift(1:3, [0, -0, 0])), ...
      permute(smoothing_filter, circshift(1:3, [0, -1, 0])));
    fx = convn(image, x_filter, 'same');
    % fy:
    y_filter = permute(x_filter, circshift(1:3, [0, -1, 0]));
    fy = convn(image, y_filter, 'same');
    % fz:
    z_filter = permute(x_filter, circshift(1:3, [0, -2, 0]));
    fz = convn(image, z_filter, 'same');
  end
  
  if use_padding
    fx = fx(pad_size + 1:end - pad_size, pad_size + 1:end - pad_size, pad_size + 1:end - pad_size);
    fy = fy(pad_size + 1:end - pad_size, pad_size + 1:end - pad_size, pad_size + 1:end - pad_size);
    fz = fz(pad_size + 1:end - pad_size, pad_size + 1:end - pad_size, pad_size + 1:end - pad_size);
  end
  
  % 3D now:
  if use_miter
    theta = corner_angles(image, miter_radius);
    theta = abs(pi - theta);
    s = max(min(1. ./ (cos(theta * .5)), 20.), 0.);
    'mean, std, min, max of corner s'
    [mean(s(:)), std(s(:))]
    [min(s(:)), max(s(:))]
    fx = fx .* s;
    fy = fy .* s;
    fz = fz .* s;
  end
  
  if size(image, 3) == 1
    %clear fz
    fz = fz .* 0;
  end
  