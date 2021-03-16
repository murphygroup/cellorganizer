function [result] = struct_ismember(a, b)
% Simple clone of ismember for structures.
% 
% Tests:
% struct_ismember([struct('a', 2), struct('a', 3)], [struct('a', 1), struct('a', 2), struct('a', 3)])
% struct_ismember([struct('a', 2), struct('a', 4)]', [struct('a', 1), struct('a', 2), struct('a', 3)])
% 
% Dependencies:
% From the File Exchange:
% cellstructeq <http://www.mathworks.com/matlabcentral/fileexchange/27542-compare-nested-cell-struct-arrays-recursively>

  % a, b, pause
  
  result = false(size(a));
  for a_index = 1:numel(a)
    for b_index = 1:length(b)
      % a(a_index), b(a_index)
      % structeq(a(a_index), b(b_index))
      if structeq(a(a_index), b(b_index))
        result(a_index) = true;
        break
      end
    end
  end

