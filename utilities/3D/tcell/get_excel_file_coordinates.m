function [result_coordinates, result_coordinate_format] = get_excel_file_coordinates(given_worksheet_data, options)
  % Returns a cell array where index is timeframe, and each element contains an array where a row is a cell (not necessarily the same one at the same index across timeframes) and columns are X, Y, T (timeframe relative to synapse formation (T=0), not evenly spaced in real time), and cell index (note: computed here, not a unique identifier if given_worksheet_data changes!).
  %
  % Dependencies:
  % File Exchange:
  % 
  % 2013-07-30 tebuck: Copied from segment_brightfield_image_by_edges.m.
  % 2014-01-22 tebuck: Check the first non-blank row to determine if this is an old-style manual one-point synapse coordinate file or a new, somewhat automatically generated (by ImageJ) two-point one. Two-point coordinates are reduced to one-point in the return value result_coordinates but not in result_two_point_coordinates.
  % 2014-03-02 tebuck: Return result_coordinates with different formats depending on annotation format instead of second output argument result_two_point_coordinates. Specify name of format in new second output argument result_coordinate_format.
  % 2014-03-24 tebuck: NOT YET: Now using indices for columns of interest (these change between old and new two-point annotations).
  
  % error
  
  default_options = struct();
  default_options.verbose = 0;
  % default_options.verbose = 1;
  

  if ~exist('options', 'var')
    options = default_options;
  else
    options = process_options_structure(default_options, options);
  end
  % options

  result_coordinates = cell(0, 1);
  file_format = '';
  cell_index = 1;
  
  % For 'new' format:
  cell_index_changed = false;
  last_line_was_header = false;
  cell_onset_time = [];
  relative_time_index = [];
  % Indices for columns of interest (these change between old and new two-point annotations):
  region_label_index = [];
  segment_distance_index = [];
  segment_angle_index = [];
  segment_left_index = [];
  segment_top_index = [];
  segment_width_index = [];
  segment_height_index = [];
  time_to_onset_index = [];
  cell_coordinate_index = [];
  onset_time_index = [];
  header_finding_function = @(given_row, given_pattern)find(~cellfun(@isempty, regexpi(given_row, given_pattern)));
  
  
  for row_index = 1:size(given_worksheet_data, 1)
    
    row = given_worksheet_data(row_index, :);
    
    % Remove spaces from cells:
    row2 = row;
    for element_index = 1:length(row2)
      if ischar(row2{element_index})
        row2{element_index} = strtrim(row2{element_index});
      end
      if isempty(row2{element_index})
        row2{element_index} = nan;
      end
    end
    row = row2;
    
    switch file_format
      case ''
        
        if ~isempty(row{1})
          if row{1}(1) == '#'
            file_format = 'old';
            result_coordinate_format = 'excel_one_point';
            % cell_index = 1;
          % elseif strcmpi(row{1}, 'Region Label')
          elseif any(strcmpi(row, 'Region Label'))
            file_format = 'new';
            result_coordinate_format = 'excel_two_point';
            last_line_was_header = true;
            % Indices for columns of interest (these change between old and new two-point annotations):
            row_string_value_indices = find(cellfun(@ischar, row));
            row_string_values = row(row_string_value_indices);
            region_label_index = row_string_value_indices(header_finding_function(row_string_values, '^Region Label$'));
            segment_distance_index = row_string_value_indices(header_finding_function(row_string_values, '^Distance$'));
            segment_angle_index = row_string_value_indices(header_finding_function(row_string_values, '^Angle$'));
            segment_left_index = row_string_value_indices(header_finding_function(row_string_values, '^Left$'));
            segment_top_index = row_string_value_indices(header_finding_function(row_string_values, '^Top$'));
            segment_width_index = row_string_value_indices(header_finding_function(row_string_values, '^Width$'));
            segment_height_index = row_string_value_indices(header_finding_function(row_string_values, '^Height$'));
            time_to_onset_index = row_string_value_indices(header_finding_function(row_string_values, '^(time to onset|relative time|realtive time)$'));
            cell_coordinate_index = row_string_value_indices(header_finding_function(row_string_values, '^(cell coordinate|coordinates|coordinate|coodinates)$'));
            onset_time_index = row_string_value_indices(header_finding_function(row_string_values, '^(onset time|onset)$'));
            if false
              % Debug info:
              row, region_label_index, segment_distance_index, segment_distance_index, segment_angle_index, segment_left_index, segment_top_index, segment_width_index, segment_height_index, time_to_onset_index, cell_coordinate_index, onset_time_index
              % keyboard
              % pause
            end
            if any(cellfun(@isempty, {region_label_index, segment_distance_index, segment_distance_index, segment_angle_index, segment_left_index, segment_top_index, segment_width_index, segment_height_index, time_to_onset_index, cell_coordinate_index, onset_time_index}))
              error('At least one required column name not found in synapse coordinates!')
            end
          else
            error('Unrecognized synapse coordinate data format!')
          end
        end
        
      case 'old'
        
        relative_time_index = row{1};
        x = row{2};
        y = row{3};
        frame_number = row{4};
        
        if ischar(relative_time_index)
          if length(relative_time_index) == 0
            % Blank line:
          elseif relative_time_index(1) == '#'
            % New cell:
            file_format = 'old';
            cell_index = cell_index + 1;
          else
            error('Line with non-numeric, non-header, and non-blank cell encountered!')
          end
          continue
        end
        
        if any(cellfun(@isnan, row(1:4)))
          % if default_options.verbose > 0, fprintf('Skipping nan line\n'), end
          continue
        end
        
        if size(result_coordinates, 1) < frame_number || isempty(result_coordinates{frame_number, 1})
          % if default_options.verbose > 0, fprintf('Initializing frame %d\n', frame_number), end
          result_coordinates{frame_number, 1} = zeros(0, 4);
        end
        
        % if default_options.verbose > 0, fprintf('Adding coordinate!\n'), end
        result_coordinates{frame_number, 1}(end + 1, :) = [x, y, relative_time_index, cell_index];
      
      case 'new'
        
        region_label = row{region_label_index};
        segment_distance = row{segment_distance_index};
        segment_angle = row{segment_angle_index};
        segment_left = row{segment_left_index};
        segment_top = row{segment_top_index};
        segment_width = row{segment_width_index};
        segment_height = row{segment_height_index};
        time_to_onset = row{time_to_onset_index};
        cell_coordinate = row{cell_coordinate_index};
        onset_time = row{onset_time_index};
        
        segment_angle = deg2rad(segment_angle);
        
        line_blank = all(cellfun(@(x)(isscalar(x) && isnan(x)), row));
        
        if line_blank || last_line_was_header
        % if isnan(region_label) || strcmpi(region_label, 'Region Label')
          % New cell, additional blank lines:
          if ~cell_index_changed
            if ~last_line_was_header
              cell_index = cell_index + 1;
            end
            cell_index_changed = true;
            cell_onset_time = [];
            relative_time_index = [];
            last_line_was_header = false;
            if default_options.verbose > 0, fprintf('cell_index = %5d\n', cell_index), end
          else
            % if default_options.verbose > 0, fprintf('Skipping blank line!\n'), end
          end
          % continue
        end
        
        % cell_index_changed = false;
        coordinates_unavailable = (isscalar(region_label) && isnan(region_label)) || any(cellfun(@(x)any(strcmpi(x, {'N/A', 'NA'})), row));
        
        if isfinite(onset_time)
          % First informative line:
          % Assume that relative_time_index, instead of being given like in 'old' files, always starts at -2 in the row where onset_time is given and goes up with every row:
          relative_time_index = -2;
          cell_onset_time = onset_time;
          if default_options.verbose > 0, fprintf('cell_onset_time = %-2d\n', cell_onset_time), end
        end
        
        if ~isempty(relative_time_index)
          if ~cell_index_changed
            % Assume that relative_time_index, instead of being given like in 'old' files, always starts at -2 in the row where onset_time is given and goes up with every row:
            relative_time_index = relative_time_index + 1;
          end
          cell_index_changed = false;
        end
        
        if ~coordinates_unavailable && ~isfinite(segment_distance)
          error('Problematic line encountered!')
        end
        
        if coordinates_unavailable
          % Otherwise uninformative line:
          if default_options.verbose > 0, fprintf('Skipping otherwise uninformative line\n'), end
          continue
        end
        
        if default_options.verbose > 0, fprintf('time_to_onset = %-2d\n', time_to_onset), end
        frame_number = cell_onset_time + time_to_onset;
        if (ischar(time_to_onset) && strcmpi(time_to_onset, 'late')) || frame_number <= 0 || isnan(time_to_onset)
          continue
        end
        
        if isnan(cell_onset_time)
          error('cell_onset_time = %-2d', cell_onset_time)
        end
        
        if default_options.verbose > 0, fprintf('frame_number = %3d\n', frame_number), end
        
        if size(result_coordinates, 1) < frame_number || isempty(result_coordinates{frame_number, 1})
          if default_options.verbose > 0, fprintf('Initializing result_coordinates{%3d, 1}\n', frame_number), end
          result_coordinates{frame_number, 1} = zeros(0, 8);
        end
        
        x = segment_left + segment_width / 2;
        y = segment_top + segment_height / 2;
        % Round so the code expecting integer coordinates works properly:
        x_unrounded = x;
        y_unrounded = y;
        x = round(x);
        y = round(y);
        if default_options.verbose > 0, fprintf('Adding coordinate x = %4d, y = %4d, cell_onset_time = %-2d, time_to_onset = %-2d, cell_index = %5d\n', x, y, cell_onset_time, time_to_onset, cell_index), end
        segment_offset = [cos(-segment_angle), sin(-segment_angle)] .* (segment_distance .* .5);
        result_coordinates{frame_number, 1}(end + 1, :) = [x, y, x_unrounded - segment_offset(1), y_unrounded - segment_offset(2), x_unrounded + segment_offset(1), y_unrounded + segment_offset(2), relative_time_index, cell_index];
        
        if false
          cell_index, cell_index_changed
          region_label, segment_distance, segment_angle, segment_left, segment_top, segment_width, segment_height, time_to_onset, cell_coordinate, onset_time, , 
          keyboard
        end
        
    end
    
    
  end
  
  % keyboard
  % error
  
  
end


