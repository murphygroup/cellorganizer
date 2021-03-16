classdef WindowedImage
  %WindowedImage Represents both compressed and uncompressed windowed images.
  % Designed to work with images of dimensionality 2 and 3 for now.

  % properties (Dependent, SetAccess = 'public')
    
  % end
  
  %grj - added functionality to specify output size by padding image on
  %WindowedImage.get_data() call

  
  properties (SetAccess = 'private')
    % Cell array with compressed windows.
    image
    % Total image size.
    image_size
    % Size of each window. Each element must be a divisor of image_size's respective element.
    window_size
    % Any one of {'none', 'gzip'}.
    compression_method
    
%     get_data_padding
  end
  
  
  methods (Access = 'public')
    function object = WindowedImage(image, window_size, compression_method)
      %WindowedImage create a windowed image either of a specified size or from an
      % existing image.
      %    If image has only one dimension with size greater than one, it is the size
      %    of a blank WindowedImage to be created. Otherwise, it is image data to be
      %    copied into a new WindowedImage.
      % get_data_padding is a 2x3 specifying the padding in the pre and
      % post directions to be applied to the image upon the get_data() call
      create_blank_image = sum(size(image) > 1) <= 1;
      if create_blank_image
        image_size = image;
      else
        image_size = size(image);
      end
      if length(image_size) < 3
        % image_size = [image_size(:), ones(1, 3 - length(image_size))]; 
        image_size = [reshape(image_size, 1, []), ones(1, 3 - length(image_size))]; 
      end
      if length(window_size) < 3
        % window_size = [window_size(:), ones(1, 3 - length(window_size))]; 
        window_size = [reshape(window_size, 1, []), ones(1, 3 - length(window_size))]; 
      end
      if ~exist('compression_method', 'var')
        compression_method = 'none'; 
      end
      if ~ismember(compression_method, {'none', 'gzip'})
        error('Compression method must be  in {''none'', ''gzip''}')
      end
      if any(mod(image_size, window_size) > 0)
        error('Image size must be divisible by window size')
      end
      
      
      object.image_size = image_size;
      object.window_size = window_size;
      object.compression_method = compression_method;
      object.image = [];
      if create_blank_image
        % % object.image = zeros(image_size([2, 1, 3:end]));
        % blank_window = object.compress_window(zeros(window_size));
        % object.image = repmat({blank_window}, image_size ./ window_size);
        object = object.set_data(0);
      else
        % keyboard
        object = object.set_data(image);
      end
      % object_post_WindowedImage = object
      % object
      % size(object.get_data()), pause
    end

    
    function window = get_array_window(a, b, window_index)
      % Gets a window from a regular array that has the same position as it
      % would in this WindowedImage.
      
      % whos
      % size(a.image)
      % window_index
      a_indices = cell(1, max(length(a.image_size), 3));
      % [a_indices{:}] = ind2sub(a.image_size, window_index);
      [a_indices{:}] = ind2sub(size(a.image), window_index);
      % a_indices
      a_indices = cell2mat(a_indices);
      lower_coordinates = a.window_size .* (a_indices - 1) + 1;
      upper_coordinates = lower_coordinates + a.window_size - 1;
      coordinates = arrayfun(...
        @(dimension_index)lower_coordinates(dimension_index):upper_coordinates(dimension_index)...
        , 1:length(lower_coordinates), 'UniformOutput', false);
      % whos coordinates
      % window = subsref(b, struct('type', '()', 'subs', coordinates));
      % Extra curly braces because of how struct treats cell arrays:
      window = subsref(b, struct('type', '()', 'subs', {coordinates}));
      % window = subsref(b, struct('type', '()', 'subs', {1, 1}));
      % a_indices
      % lower_coordinates
      % upper_coordinates
      % coordinates
      % whos window
    end
    
    
    function window = get_window(a, b, window_index, pad_radius, mode)
      % Gets a window from this WindowedImage with optional padding (which
      % optionally treats space as periodic).

      error('This more generic version would be useful, but use get_padded_window for now')
      
      if length(pad_radius) < length(a.image_size)
        pad_radius = [pad_radius, zeros(1, length(a.image_size) - length(pad_radius))]; 
      end
      
      dimensions_to_pad = find(pad_radius > 0); 
      number_windows_to_pad = ceil(pad_radius ./ a.window_size); 
      
      % whos
      % size(a.image)
      % window_index
      a_indices = cell(1, max(length(a.image_size), 3));
      % [a_indices{:}] = ind2sub(a.image_size, window_index);
      [a_indices{:}] = ind2sub(size(a.image), window_index);
      % a_indices
      a_indices = cell2mat(a_indices);
      lower_coordinates = a.window_size .* (a_indices - 1) + 1;
      upper_coordinates = lower_coordinates + a.window_size - 1;
      coordinates = arrayfun(...
        @(dimension_index)lower_coordinates(dimension_index):upper_coordinates(dimension_index)...
        , 1:length(lower_coordinates), 'UniformOutput', false);
      % whos coordinates
      % window = subsref(b, struct('type', '()', 'subs', coordinates));
      % Extra curly braces because of how struct treats cell arrays:
      window = subsref(b, struct('type', '()', 'subs', {coordinates}));
      % window = subsref(b, struct('type', '()', 'subs', {1, 1}));
      % a_indices
      % lower_coordinates
      % upper_coordinates
      % coordinates
      % whos window
    end
    
    
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
      
      [number_rows, number_columns] = size(a.image); 
      %whos
      %row, column
      padded_window = padarray(...
        a.decompress_window(a.image{row, column}),...
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
        adjacent_image = a.decompress_window(...
          a.image{adjacent_row, adjacent_column});
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
    
    
    function c = set_window(a, row, column, b)
      c = a;
      c.image{row, column} = c.compress_window(b);
    end
    
    
    function image = get_data(a)
      % a
      % a.image
      image = cell2mat(cellfun(@(index)a.decompress_window(index), a.image, 'UniformOutput', false));
      
%       final_pad = a.get_data_padding;
%       
%       final_pad_pos = final_pad;
%       final_pad_neg = -final_pad;
%       
%       final_pad_pos(final_pad < 0) = 0;
%       final_pad_neg(final_pad > 0) = 0;
%       
%       image = padarray(image, final_pad_pos(1,:), 0, 'pre');
%       image = padarray(image, final_pad_pos(2,:), 0, 'post');
%       
%       image = image((1+final_pad_neg(1,1)):(end-final_pad_neg(2,1)), (1+final_pad_neg(1,2)):(end-final_pad_neg(2,2)),(1+final_pad_neg(1,3)):(end-final_pad_neg(2,3)));
%       
    end
    
    
    function c = set_data(a, b)
      % if ~(((numel(a.image_size) == numel(size(b))) && all(a.image_size == size(b))) || (numel(b) == 1))
      if ~(((sum(a.image_size > 1) == numel(size(b))) && all(a.image_size(a.image_size > 1) == size(b))) || (numel(b) == 1))
        error('Numeric operand must be scalar or the same size as the image')
      end
      % if ~((strcmp(class(a), class(b)) && numel(a.image_size) == numel(size(b)) && all(a.image_size == size(b))) || (isnumeric(b) && numel(b) == 1))
        % error('Numeric operand must be scalar or the same size as the image')
      % end
      % a
      % whos
      if numel(b) == 1
        blank_window = a.compress_window(repmat(b, a.window_size));
        a.image = repmat({blank_window}, a.image_size ./ a.window_size);
        % a.image
      else
        a.image = cell(a.image_size ./ a.window_size);
        % a.image
        % error
        for window_index = 1:numel(a.image)
          b_window = a.get_array_window(b, window_index);
          % whos b_window, pause
          a.image{window_index} = a.compress_window(b_window);
        end
        % a.image
      end
      % a_post_set_data = a
      c = a;
      % error
    end
    
    
    function c = set_borders(a, b, only_xy)
      % Only works with 2D images or 3D images with window_size(3) == image_size(3).
      if ~exist('only_xy', 'var')
        only_xy = false; 
      end

      if ~(isnumeric(b) && numel(b) == 1)
        error('Numeric operand must be scalar')
      end
      
      [number_rows, number_columns] = size(a.image);
      for row = 1:number_rows
        for column = 1:number_columns
          if all(row ~= [1, number_rows]) && all(column ~= [1, number_columns])
            continue
          end
          
          window = a.decompress_window(a.image{row, column});
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
          
          a.image{row, column} = a.compress_window(window);
        end
      end
      
      c = a;
      
      % for window_index = 1:numel(c.image)
        % b_window = c.get_array_window(b, window_index);
        % c.image{window_index} = c.compress_window(b_window);
      % end
    end
    
    
    function c = binary_operation(a, b, function_handle)
      % binary_operation(a, b, function_handle)
      % Helper function called by window-wise (usually element-wise) functions like plus, minus, and times.
      
      % % function object = binary_operation(a, b, function_handle, should_be_same_size)
      % if ~exist('should_be_same_size', 'var')
        % should_be_same_size = true;
      % end
      
      % If one is numeric, assume the other is a WindowedImage:
      if (isnumeric(a) || islogical(a)) && numel(a) > 1
        a = WindowedImage(a, b.window_size, b.compression_method);
      elseif (isnumeric(b) || islogical(b)) && numel(b) > 1
        b = WindowedImage(b, a.window_size, a.compression_method);
      end
      
      if strcmp(class(a), class(b))
        if a.image_size ~= b.image_size
          error('Images must be the same size')
        end
        if a.window_size ~= b.window_size
          error('Images'' windows must be the same size')
        end
      elseif (isnumeric(b) || islogical(b))
      else
        % whos
        error('Operands that are not images must be numeric')
      end
      % c = WindowedImage(a.image_size, a.window_size, a.compression_method);
      c = a;
      
      
      for window_index = 1:numel(a.image)
        if strcmp(class(a), class(b))
          c.image{window_index} = c.compress_window(function_handle(...
            a.decompress_window(a.image{window_index}), b.decompress_window(b.image{window_index}) ...
            ));
        elseif isnumeric(b) && numel(b) > 1
          % a_indices = cell(1, ndims(a.image));
          % [a_indices{:}] = ind2sub(size(a.image), window_index);
          % a_indices = cell2mat(a_indices);
          % lower_coordinates = a.window_size .* (a_indices - 1) + 1;
          % upper_coordinates = lower_coordinates + a.window_size - 1;
          % coordinates = arrayfun(...
            % @(dimension_index2)lower_coordinates(dimension_index2):upper_coordinates(dimension_index2)...
            % , 1:ndims(lower_coordinates), 'UniformOutput', false);
          % b_window = subsref(b, struct('type', '()', 'subs', coordinates));
          b_window = a.get_array_window(b, window_index);
          c.image{window_index} = c.compress_window(function_handle(...
            a.decompress_window(a.image{window_index}), b_window ...
            ));
        else
          c.image{window_index} = c.compress_window(function_handle(...
            a.decompress_window(a.image{window_index}), b ...
            ));
        end
      end
    end

    % % Would have been cool if this had worked:
    % function_list = {@plus, @minus, @times, @rdivide, @ldivide, @power, @lt, @gt, @le, @ge, @ne};
    % for f = function_list
      % fs = func2str(f);
      % eval(fprintf('function c = %s(a, b), c = a.binary_operation(b, @%s); end', fs, fs))
    % end
    
    function c = plus(a, b), c = binary_operation(a, b, @plus); end
    function c = minus(a, b), c = binary_operation(a, b, @minus); end
    function c = times(a, b), c = binary_operation(a, b, @times); end
    function c = rdivide(a, b), c = binary_operation(a, b, @rdivide); end
    function c = ldivide(a, b), c = binary_operation(a, b, @ldivide); end
    function c = power(a, b), c = binary_operation(a, b, @power); end
    function c = eq(a, b), c = binary_operation(a, b, @eq); end
    function c = ne(a, b), c = binary_operation(a, b, @ne); end
    function c = lt(a, b), c = binary_operation(a, b, @lt); end
    function c = gt(a, b), c = binary_operation(a, b, @gt); end
    function c = le(a, b), c = binary_operation(a, b, @le); end
    function c = ge(a, b), c = binary_operation(a, b, @ge); end
    function c = and(a, b), c = binary_operation(a, b, @and); end
    function c = or(a, b), c = binary_operation(a, b, @or); end
    function c = xor(a, b), c = binary_operation(a, b, @xor); end
    % function c = max(a, b), c = binary_operation(a, b, @max); end
    % function c = min(a, b), c = binary_operation(a, b, @min); end
    
    % function c = mtimes(a, b), error('Not implemented'); end
    % Also need:
    % abs, max, min, set_borders, size, sum, negate?


    
    function b = unary_operation(a, function_handle)
      % unary_operation(a, function_handle)
      % Helper function called by window-wise (usually element-wise) functions like uminus and abs.

      b = WindowedImage(a.image_size, a.window_size, a.compression_method);
      
      for window_index = 1:numel(a.image)
        b.image{window_index} = b.compress_window(function_handle(...
          a.decompress_window(a.image{window_index})...
          ));
      end
    end

    function b = uplus(a), b = a.unary_operation(@uplus); end
    function b = uminus(a), b = a.unary_operation(@uminus); end
    function b = not(a), b = a.unary_operation(@not); end
    % function b = any(a), b = a.unary_operation(@any); end
    % function b = all(a), b = a.unary_operation(@all); end
    function b = abs(a), b = a.unary_operation(@abs); end
    
    
    function b = summarizing_unary_operation(a, function_handle)
      % unary_operation(a, function_handle)
      % Helper function called by element-wise summarizing functions like sum and min.
      b = unary_operation(a, @(value)function_handle(value(:)));
      b = cellfun(@(value)b.decompress_window(value), b.image);
      b = function_handle(b(:));
    end

    
    % Does not support operating along a dimension yet.
    function b = sum(a), b = a.summarizing_unary_operation(@sum); end
    function b = any(a), b = a.summarizing_unary_operation(@any); end
    function b = all(a), b = a.summarizing_unary_operation(@all); end

    % function b = min(a), b = a.summarizing_unary_operation(@min); end
    % function b = max(a), b = a.summarizing_unary_operation(@max); end
    

    function c = min(a, b)
      % Does not support operating along a dimension yet.
      h = @min;
      if exist('b', 'var')
        c = binary_operation(a, b, h);
      else
        c = a.summarizing_unary_operation(h);
      end
    end

    function c = max(a, b)
      % Does not support operating along a dimension yet.
      h = @max;
      if exist('b', 'var')
        c = binary_operation(a, b, h);
      else
        c = a.summarizing_unary_operation(h);
      end
    end

    
    
    % Nevermind, just provide the necessary functions in this class instead of allowing indexing:
    % function varargout = subsref(object, index_structure)
      % if length(index_structure) > 1
        % error('Multiple indexing expressions not yet implemented')
      % end
      % if strcmp(index_structure.type, '.')
        % % From http://stackoverflow.com/questions/8852603/access-to-methods-and-properties-for-classes-with-overloaded-index-reference
        % varargout = cell(1, nargout);
        % if nargout == 0
          % builtin('subsref', object, index_structure);
        % else
          % [varargout{:}] = builtin('subsref', object, index_structure);
        % end
      % end
      % for dimension_index = 1:ndims(object.image)
        % if strcmp(index_structure.subs{dimension_index}, ':')
      % result = cellfun()
    % end

    
    function c = interp3(a, b, method, extrapval, pad_radius, window_mode)
      % interp3(a, b, method, extrapval, pad_radius, window_mode)
      % Finds values of a at locations b like Matlab's interp3.

      if iscell(b) && numel(b) == 3 && all(cellfun(@(value)strcmp(class(a), class(value)), b))
        if any(cellfun(@(value)any(a.image_size ~= value.image_size), b))
          error('Images must be the same size')
        end
        if any(cellfun(@(value)any(a.window_size ~= value.window_size), b))
          error('Images'' windows must be the same size')
        end
      elseif isnumeric(b)
        if ~((...
          all(cellfun(@(value)numel(a.image_size) == numel(size(value)), b)) && ...
          all(cellfun(@(value)all(a.image_size == size(value)), b))) || ...
          all(cellfun(@(value)numel(value) == 1, b)))
          error('Numeric operand must be scalar or the same size as the image')
        end
      else
        error('Operands that are not images must be numeric')
      end

      if ~exist('method', 'var')
        method = '*linear';
      end
      if ~exist('extrapval', 'var')
        extrapval = nan;
      end
      if ~exist('window_mode', 'var')
        window_mode = 'circular';
      end
      
      decompress_whole = ~exist('pad_radius', 'var'); 
      if ~decompress_whole && pad_radius < 0
        decompress_whole = true;
      end
      
      % c = WindowedImage(a.image_size, a.window_size, a.compression_method);
      c = a;
      
      if (decompress_whole)
        a2 = a.get_data();
        % whos
        % a
        % a.image
        for ind = 1:numel(a.image)
          if ndims(a2) == 3
            c.image{ind} = c.compress_window(interp3(...
              a2, ...
              b{1}.decompress_window(b{1}.image{ind}), ...
              b{2}.decompress_window(b{2}.image{ind}), ...
              b{3}.decompress_window(b{3}.image{ind}), ...
              method, extrapval));
          elseif ndims(a2) == 2
            c.image{ind} = c.compress_window(interp2(...
              a2, ...
              b{1}.decompress_window(b{1}.image{ind}), ...
              b{2}.decompress_window(b{2}.image{ind}), ...
              method, extrapval));
          end
          
        end
      else
        error('Not yet implemented')
        % a2 = a;
        % [number_rows, number_columns] = size(a2); 
        % %b{1, 1}
        % window_size = size(CompressLib.decompress(b{1}{1, 1}));
        % for row = 1:number_rows
          % for column = 1:number_columns
            % window = get_compressed_window(...
              % a2, row, column, pad_radius, window_mode); 
            % if numel(window_size) == 3
              % a{row, column} = CompressLib.compress(interp3(...
                % window, ...
                % CompressLib.decompress(b{1}{row, column})...
                % - window_size(2) * (column - 1) + pad_radius, ...
                % CompressLib.decompress(b{2}{row, column})...
                % - window_size(1) * (row - 1) + pad_radius, ...
                % CompressLib.decompress(b{3}{row, column}), ...
                % method, extrapval));
            % else
              % a{row, column} = CompressLib.compress(interp2(...
                % window, ...
                % CompressLib.decompress(b{1}{row, column})...
                % - window_size(2) * (column - 1) + pad_radius, ...
                % CompressLib.decompress(b{2}{row, column})...
                % - window_size(1) * (row - 1) + pad_radius, ...
                % method, extrapval));
            % end
          % end
        % end
      end
    end
  
  
  end
  
  
  methods (Access = 'private')
    function compressed_window = compress_window(object, individual_window)
      compressed_window = individual_window;
      if strcmp(object.compression_method, 'gzip')
        compressed_window = CompressLib.compress(compressed_window);
      end
    end

    function decompressed_window = decompress_window(object, individual_window)
      decompressed_window = individual_window;
      if strcmp(object.compression_method, 'gzip')
        decompressed_window = CompressLib.decompress(decompressed_window);
      end
    end
  end
  
  
  methods (Static, Access = 'public')
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
      [number_rows, number_columns] = size(a{1}.image); 
      %b{1, 1}
      window_size = a{1}.window_size;
      for row = 1:number_rows
        for column = 1:number_columns
          Vx1 = a{1}.get_padded_window(row, column, 1, window_mode); 
          Vy1 = a{2}.get_padded_window(row, column, 1, window_mode); 
          Vz1 = a{3}.get_padded_window(row, column, 1, window_mode); 
          dx2 = (-del2(Vx1)*6*alpha + gamma*Vx1).^2;
          dy2 = (-del2(Vy1)*6*alpha + gamma*Vy1).^2;
          dz2 = (-del2(Vz1)*6*alpha + gamma*Vz1).^2;
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
  
  
  

    %
    % Testing suite:
    %
    
    function test()
      % Runs all methods prefixed with "test_" (methods are assumed to be static).
      fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
      fprintf('Running tests\n')
      
      test_method_names = methods('WindowedImage');
      test_method_names = test_method_names(strncmp('test_', test_method_names, 5));
      methods_tested = {};
      for test_method_name = test_method_names(:)'
        test_method_name = test_method_name{1};
        % test_method_name
        fprintf(['Running ', test_method_name, '...\n'])
        novel_methods_tested = eval([mfilename('class'), '.', test_method_name]);
        methods_tested = [methods_tested; novel_methods_tested(:)];
      end

      methods_untested = methods('WindowedImage');
      methods_untested = methods_untested(~strncmp('test_', methods_untested, 5));
      % This method:
      methods_untested = methods_untested(~strncmp('test', methods_untested, 4));
      methods_untested = methods_untested(~ismember(methods_untested, methods_tested));
      fprintf('Untested methods:\n')
      for untested_method_name = methods_untested(:)'
        untested_method_name = untested_method_name{1};
        fprintf([' ', untested_method_name, '\n'])
      end
    end
  
  
    function methods_tested = test_binary_operations()
      methods = {@plus, @minus, @times, @rdivide, @ldivide, @power, @eq, @ne, @lt, @gt, @le, @ge, @and, @or, @xor, @min, @max};
      % which plus
      n = 100;
      p = 10;
      wn = 50;
      a = rand(n, n, p) + 1;
      b = rand(n, n, p) + 1;
      for method = methods
        method = method{1};
        % fprintf([' Testing ', func2str(method), '\n']);
        fprintf([' ', func2str(method), '\n']);
        c = method(a, b);
        g = method(WindowedImage(a, [wn, wn, p], 'none'), WindowedImage(b, [wn, wn, p], 'none'));
        h = method(WindowedImage(a, [wn, wn, p], 'gzip'), WindowedImage(b, [wn, wn, p], 'gzip'));
        condition1 = get_data(c ~= g);
        condition1 = any(condition1(:));
        condition2 = get_data(c ~= h);
        condition2 = any(condition2(:));
        if condition1 || condition2
        % if condition1
          % error(['WindowedImage implementation for ', func2str(method), ' incorrect!'])
          error(['WindowedImage implementation incorrect!'])
        end
      end
      
      methods_tested = cellfun(@(x)func2str(x), methods, 'UniformOutput', false);
      methods_tested{end+1} = 'binary_operation';
    end

    function methods_tested = test_unary_operations()
      methods = {@uplus, @uminus, @not, @abs};
      % which plus
      n = 100;
      p = 10;
      wn = 50;
      a = rand(n, n, p) + 1;
      for method = methods
        method = method{1};
        fprintf([' ', func2str(method), '\n']);
        c = method(a);
        g = method(WindowedImage(a, [wn, wn, p], 'none'));
        h = method(WindowedImage(a, [wn, wn, p], 'gzip'));
        condition1 = get_data(c ~= g);
        condition1 = any(condition1(:));
        condition2 = get_data(c ~= h);
        condition2 = any(condition2(:));
        if condition1 || condition2
          error(['WindowedImage implementation incorrect!'])
        end
      end
      
      methods_tested = cellfun(@(x)func2str(x), methods, 'UniformOutput', false);
      methods_tested{end+1} = 'unary_operation';
    end
    
    function methods_tested = test_summarizing_unary_operations()
      methods = {@sum, @any, @all, @min, @max};
      % which plus
      n = 100;
      p = 10;
      wn = 50;
      a = rand(n, n, p) + 1;
      for method = methods
        method = method{1};
        % fprintf([' Testing ', func2str(method), '\n']);
        fprintf([' ', func2str(method), '\n']);
        c = method(a(:));
        g = method(WindowedImage(a, [wn, wn, p], 'none'));
        h = method(WindowedImage(a, [wn, wn, p], 'gzip'));
        % whos
        % condition1 = c ~= g;
        % condition2 = c ~= h;
        condition1 = abs(c - g) >= 1e-6;
        condition2 = abs(c - h) >= 1e-6;
        % c, g, h
        % abs(c - g)
        % abs(c - h)
        if condition1 || condition2
        % if condition1
          % error(['WindowedImage implementation for ', func2str(method), ' incorrect!'])
          error(['WindowedImage implementation incorrect!'])
        end
      end
      
      methods_tested = cellfun(@(x)func2str(x), methods, 'UniformOutput', false);
      methods_tested{end+1} = 'summarizing_unary_operation';
    end
    
  end
  
end