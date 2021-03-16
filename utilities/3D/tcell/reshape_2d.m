function [a] = reshape_2d(a, r)
  % Reshape a 3D image a so that its Z slices (along the third dimension) become arranged left to right in a 2D image. If r > 1, wrap the slices like left-to-right text into r rows. If r == -1, automatically determine the number of slices to get an approximately square image.
  if ~exist('r', 'var') || isempty(r)
    r = 1;
  elseif r == -1
    r = max(floor(sqrt(size(a, 3) * size(a, 2) / size(a, 1))), 1);
  end
  s = size(a);
  c = round(s(3) / r);
  b = zeros(s(1) * r, s(2) * c) + mean(a(:));
  for row_index = 1:r
    for column_index = 1:c
      if (row_index - 1) * c + column_index > s(3)
        break
      end
      rs = s(1) * (row_index - 1) + 1;
      rf = s(1) * (row_index);
      cs = s(2) * (column_index - 1) + 1;
      cf = s(2) * (column_index);
      b(rs:rf, cs:cf) = a(:, :, (row_index - 1) * c + column_index);
    end
  end
  a = b;
  
