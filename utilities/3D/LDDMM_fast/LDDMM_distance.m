    function d = LDDMM_distance(a, alpha, gamma, window_mode)
      % Input is a cell array of three compressed images representing
      % the X, Y, and Z velocities of an LDDMM integration step.
      % Assumes input has zero borders with window_mode = 'replicate'.
      if ~exist('window_mode', 'var')
        %window_mode = 'circular';
        window_mode = 'replicate';
      end
      if ~exist('alpha', 'var')
        alpha = 1.;
      end
      if ~exist('gamma', 'var')
        gamma = .04;
      end

      d = 0;
      [number_rows, number_columns] = size(a{1}); 
      %b{1, 1}
     %  window_size = a{1}.window_size;
      for row = 1:number_rows
        for column = 1:number_columns
          Vx1 = get_padded_window(a{1}, row, column, 1, window_mode); 
          Vy1 = get_padded_window(a{2}, row, column, 1, window_mode); 
          Vz1 = get_padded_window(a{3}, row, column, 1, window_mode); 
          % xruan 09/06/2015
          %dx2 = (-del2(Vx1)*6*alpha + gamma*Vx1).^2;
          %dy2 = (-del2(Vy1)*6*alpha + gamma*Vy1).^2;
          %dz2 = (-del2(Vz1)*6*alpha + gamma*Vz1).^2;
          dx2 = (-del2_2d(Vx1)*6*alpha + gamma*Vx1).^2;
          dy2 = (-del2_2d(Vy1)*6*alpha + gamma*Vy1).^2;
          dz2 = (-del2_2d(Vz1)*6*alpha + gamma*Vz1).^2;
          dx2 = dx2(2:end-1, 2:end-1, :);
          dx2 = mean(dx2(:));
          dy2 = dy2(2:end-1, 2:end-1, :);
          dy2 = mean(dy2(:));
          dz2 = dz2(2:end-1, 2:end-1, :);
          dz2 = mean(dz2(:));
          %d = d + sqrt(dx2 + dy2 + dz2);
          d = d + dx2 + dy2 + dz2;
        end
      end
      % Why is this scaled by the number of windows?????
      d = sqrt(d / (number_rows * number_columns * 1.)); 
    end