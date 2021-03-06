    function padded_window = get_padded_window(a, row, column, pad_radius, mode)
      % For each window, pad it with paddarray in replicate mode with
      % pad_radius + 1 elements in the first two dimensions and
      % replace the padding with portions of the adjacent windows, if
      % possible. Run single_step_window on the padded window, store
      % its returned derivatives in compressed form, and compute the
      % distance moved by the deformation. Note that this assumes that
      % windows are larger than pad_radius.
      
      if (nargin < 5)
        %mode = 'replicate';
        mode = 'circular';
      end
      
      number_rows = 1;
       number_columns = 1;
      %whos
      %row, column
      padded_window = padarray(...
        a,...
        [1, 1, 0] * (pad_radius), ...
        'both', 'replicate'); 
      if (pad_radius == 0)
        return 
      end
      
      % Offsets for rows and columns, starting from the north side and
      % going clockwise to the northwest side:
      edge_regions = [...
        -1, 0; ...
        -1, 1; ...
        0, 1; ...
        1, 1; ...
        1, 0; ...
        1, -1; ...
        0, -1; ...
        -1, -1; ...
                     ]; 
      
      for edge_index = 1:size(edge_regions, 1)
        r = edge_regions(edge_index, 1); 
        c = edge_regions(edge_index, 2); 
        if ~(strcmp(mode, 'circular') || ...
             ((row + r > 0) && (row + r <= number_rows) && ...
              (column + c > 0) && (column + c <= number_columns)))
          continue
        end
        adjacent_row = mod(row - 1 + r, number_rows) + 1; 
        adjacent_column = mod(column - 1 + c, number_columns) + 1; 
        adjacent_image = a{adjacent_row, adjacent_column};
        % Indices for adjacent_image:
        adjacent_row_start = ...
            (r == -1) * (size(adjacent_image, 1) - pad_radius + 1) + ...
            (r == 0) * (1) + ...
            (r == 1) * (1) ...
            ;
        adjacent_row_finish = ...
            (r == -1) * (size(adjacent_image, 1)) + ...
            (r == 0) * (size(adjacent_image, 1)) + ...
            (r == 1) * (pad_radius) ...
            ;
        adjacent_column_start = ...
            (c == -1) * (size(adjacent_image, 2) - pad_radius + 1) + ...
            (c == 0) * (1) + ...
            (c == 1) * (1) ...
            ;
        adjacent_column_finish = ...
            (c == -1) * (size(adjacent_image, 2)) + ...
            (c == 0) * (size(adjacent_image, 2)) + ...
            (c == 1) * (pad_radius) ...
            ;
        % Indices for padded_window:
        window_row_start = ...
            (r == -1) * (1) + ...
            (r == 0) * (pad_radius + 1) + ...
            (r == 1) * (size(padded_window, 1) - pad_radius + 1) ...
            ;
        window_row_finish = ...
            (r == -1) * (pad_radius) + ...
            (r == 0) * (size(padded_window, 1) - pad_radius) + ...
            (r == 1) * (size(padded_window, 1)) ...
            ;
        window_column_start = ...
            (c == -1) * (1) + ...
            (c == 0) * (pad_radius + 1) + ...
            (c == 1) * (size(padded_window, 2) - pad_radius + 1) ...
            ;
        window_column_finish = ...
            (c == -1) * (pad_radius) + ...
            (c == 0) * (size(padded_window, 2) - pad_radius) + ...
            (c == 1) * (size(padded_window, 2)) ...
            ;
        padded_window(...
          window_row_start:window_row_finish, ...
          window_column_start:window_column_finish, ...
          :) = adjacent_image(...
            adjacent_row_start:adjacent_row_finish, ...
            adjacent_column_start:adjacent_column_finish, ...
            :); 
      end
    end