   function c = set_borders(a, b, only_xy)
      % Only works with 2D images or 3D images with window_size(3) == image_size(3).
      if ~exist('only_xy', 'var')
        only_xy = false; 
      end

      if ~(isnumeric(b) && numel(b) == 1)
        error('Numeric operand must be scalar')
      end
      
      [number_rows, number_columns] = size(1);
      for row = 1:number_rows
        for column = 1:number_columns
          if all(row ~= [1, number_rows]) && all(column ~= [1, number_columns])
            continue
          end
          
          window = a;
          if (row == 1)
            window(1, :, :)=b; 
          end
          if (row == number_rows)
            window(end, :, :)=b; 
          end
          if (column == 1)
            window(:, 1, :)=b; 
          end
          if (column == number_columns)
            window(:, end, :)=b; 
          end
          if ~only_xy
            window(:, :, 1)=b; 
            window(:, :, end)=b; 
          end
          
          a =window;
        end
      end
      
      c = a;
      
      % for window_index = 1:numel(c.image)
        % b_window = c.get_array_window(b, window_index);
        % c.image{window_index} = c.compress_window(b_window);
      % end
    end