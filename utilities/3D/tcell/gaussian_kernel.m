function [f] = gaussian_kernel(sigma)
  % 2013-04-01 tebuck: now all dimensions of output are odd-sized.
  if numel(sigma) == 1
    sigma = ones(1, 3) * sigma;
  end

  s = ceil(sigma * 3.) * 2 + 1;
  % [cx, cy, cz] = meshgrid(1:s, 1:s, 1:s);
  [cx, cy, cz] = meshgrid(1:s(1), 1:s(2), 1:s(3));
  m = (s + 1) / 2;
  f = zeros(size(cx));
  % f = reshape(mvnpdf([cx(:), cy(:), cz(:)], ceil(size(cx)/2), diag(ones(1, 3) * sigma).^2), size(cx));
  % f = reshape(mvnpdf([cx(:), cy(:), cz(:)], ceil(size(cx)/2), diag(sigma).^2), size(cx));
  f = reshape(mvnpdf([cx(:), cy(:), cz(:)], m, sigma.^2), size(cx));
  f = f ./ sum(f(:));
  
  % % Debugging:
  % a = zeros(size(cx));
  % a(m(2) + (-1:1), m(1) + (-1:1), m(3) + (-1:1)) = 1;
  % imshow(reshape_contrast(a)), pause
  % imshow(reshape_contrast(f)), pause
