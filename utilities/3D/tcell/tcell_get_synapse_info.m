function [t_cell_info] = tcell_get_synapse_info(t_cell_info, options)
  % 2016-02-25 xruan: Copied from master_script_get_synapse_info.m.
  % 
  % 
  % 
  % 
  % % Dependencies:
  % % From the File Exchange: export_fig, pmkmp, toolbox_graph, toolbox_fast_marching, affine, Snake3D, Mesh_voxelisation, dirr, inhull
  % % From tebuck: modified version of Jieyue Li's HPA_lib, rasterize_mesh

  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Produce lists of synapse files, images, and individual synapses, frames, and cells:
  

  % Variables referenced more than once can be copied to local variables:
  master_script_options = t_cell_info.options;
  synapse_location = master_script_options.synapse_location;
  image_location = master_script_options.image_location;
  
  % sensors_to_exclude = master_script_options.sensors_to_exclude;
  % sensors_to_include = master_script_options.sensors_to_include;
  timepoints_to_include = master_script_options.timepoints_to_include;
  % conditions_to_include = master_script_options.conditions_to_include;
  
  t_cell_info.synapse_info = struct();

  
  function [result_path] = clean_relative_path(given_path)
    % Remove any initial slash and reduce repeated slashes:
    result_path = given_path;
    % result_path = regexprep(result_path, '^/', '');
    % result_path = regexprep(result_path, '/+', '/');
    % result_path = regexprep([result_path, '/'], '/+$', '/');
    % result_path = regexprep([result_path, '/'], '/+', '/');
    result_path = regexprep(strcat(result_path, '/'), '/+', '/');
    result_path = regexprep(result_path, '^/', '');
  end
  
  
  function [result_path] = clean_absolute_path(given_path)
    % Remove any initial slash and reduce repeated slashes:
    result_path = given_path;
    % result_path = regexprep([result_path, '/'], '/+', '/');
    result_path = regexprep(strcat(result_path, '/'), '/+', '/');
  end
  
  
  function [result_protein_synonyms] = get_protein_synonyms(given_synonym)
    % Return the given protein name unless master_script_options.protein_name_synonyms specifies the proper name for it:
    result_protein_synonyms = {given_synonym};
    for protein_synonym_index = 1:length(master_script_options.protein_name_synonyms)
      if any(strcmp(master_script_options.protein_name_synonyms{protein_synonym_index}, given_synonym))
        result_protein_synonyms = master_script_options.protein_name_synonyms{protein_synonym_index};
        return
      end
    end
  end
  
  
  function [result_protein_name] = get_protein_name_given_synonym(given_synonym)
    % Return the given protein name unless master_script_options.protein_name_synonyms specifies the proper name for it:
    result_protein_name = get_protein_synonyms(given_synonym);
    result_protein_name = result_protein_name{1};
  end
  
  
  t_cell_info.synapse_info.clean_relative_path = @clean_relative_path;
  t_cell_info.synapse_info.clean_absolute_path = @clean_absolute_path;
  t_cell_info.synapse_info.get_protein_synonyms = @get_protein_synonyms;
  t_cell_info.synapse_info.get_protein_name_given_synonym = @get_protein_name_given_synonym;
  
  
  % 2014-01-21 tebuck: modifying so that synapse_location can now be a cell array of locations:
  
  if ~iscell(synapse_location)
    synapse_location = {synapse_location};
  end
  
  if ~iscell(image_location)
    image_location = {image_location};
  end
  
  
  % Get a list of all synapse files:
  % Read CSV files to get synapse information, and use their filenames to search for the correct image directory and images:
  
  coordinate_file_list = [];
  for synapse_file_index = 1:length(synapse_location)
    current_synapse_location = synapse_location{synapse_file_index};
    [current_synapse_location, name, ext]= fileparts(current_synapse_location);
    current_coordinate_file_list = dir([current_synapse_location, '/', name, '.csv']);
    coordinate_file_list = [coordinate_file_list; strcat(current_synapse_location, '/', {current_coordinate_file_list.name}.')];
  end
  run_list = coordinate_file_list;
%   % only use included condition and sensor
%   run_list = {};
%   for condition_index = 1 : numel(conditions_to_include)
%         for sensor_index = 1 : numel(sensors_to_include)
%             curr_condition = conditions_to_include{condition_index};
%             curr_sensor = sensors_to_include{sensor_index};
%             condition_sensor_combination = [curr_condition, '_', curr_sensor];
%             coordinate_file_inds = ~cellfun(@isempty, regexp(coordinate_file_list, condition_sensor_combination));
%             if any(coordinate_file_inds)
%                 run_list(end + 1 : end + sum(coordinate_file_inds), 1) = coordinate_file_list(coordinate_file_inds);
%             else
%                 warning(sprintf('Condition %s - Sensor %s combination does not exist!', curr_condition, curr_sensor));
%             end  
%         end
%   end
  
  number_runs = length(run_list);

  
  % warning('>>>> HACK, setting number_runs = 1 for testing purposes'), number_runs = 1;
  
  run_relative_path_list = cell(number_runs, 1);
  run_key_list = cell(number_runs, 1);
  
  run_data_list = cell(number_runs, 1);
  run_data_format_list = cell(number_runs, 1);
  run_images_list = cell(number_runs, 1);
  run_brightfield_images_list = cell(number_runs, 1);
  run_image_frame_indices_list = cell(number_runs, 1);
  run_image_size_list = zeros(number_runs, 3);
  run_image_voxel_size_list = zeros(number_runs, 3);
  
  run_number_cells = zeros(number_runs, 1);
  run_number_frames = zeros(number_runs, 1);
  
  run_point_tracks = cell(number_runs, 1);
  run_relative_time_tracks = cell(number_runs, 1);
  
  run_relative_times = cell(number_runs, 1);
  run_number_relative_times = zeros(number_runs, 1);
  all_relative_times = [];
  number_all_relative_times = 0;
  all_relative_times_seconds_data = [];
  all_relative_times_seconds = [];
  
  % runs_to_keep = true(number_runs, 1);
  runs_to_keep = false(number_runs, 1);
  
  filterable_run_variables = {'run_list', 'run_relative_path_list', 'run_key_list', 'run_data_list', 'run_data_format_list', 'run_images_list', 'run_brightfield_images_list', 'run_image_frame_indices_list', 'run_image_size_list', 'run_image_voxel_size_list', 'run_number_cells', 'run_number_frames', 'run_point_tracks', 'run_relative_time_tracks', 'run_relative_times', 'run_number_relative_times'};
  
  
  function [result] = conditional_str2num(given_string)
    % Convert to numeric value if it is numeric, otherwise keep the string value:
    result = given_string;
    if ~isnumeric(result)
      result = str2num(result);
    end
    if isempty(result)
      result = given_string;
    end
  end
  
  function [result] = conditional_str2double(given_string)
    % Convert to numeric value if it is numeric, otherwise keep the string value (but using str2double):
    if isempty(given_string)
      result = nan;
    elseif ~isnumeric(given_string)
      result = str2double(given_string);
      if isnan(result) && ~strcmpi(given_string, 'nan')
        result = given_string;
      end
    else
      result = given_string;
    end
  end
  
  
  % Collect synapse info:
  % fprintf('\nCollecting synapse locations from Excel files in "%s"\n', master_script_options.annotation_location)
  % fprintf('\nCollecting synapse locations from Excel files in "%s"\n', master_script_options.synapse_location)
  fprintf('\nCollecting synapse locations from Excel files in the following locations:\n');
  for synapse_file_index = 1:length(synapse_location)
    fprintf('  %s\n', synapse_location{synapse_file_index})
  end
  for run_index = 1:number_runs
    % Filename of file specifying synapse locations in 2D and time:
    run_file = run_list{run_index};

    % fprintf('Collect synapse info: %s\n', run_file)
    
    [~, run_file_name, run_file_extension] = fileparts(run_file);
    run_file_name_tokens = regexp(run_file_name, '(.+)_(.+)_(.+)', 'tokens');
    run_base_subdirectory = run_file_name_tokens{1}{1};
    run_protein = run_file_name_tokens{1}{2};
    run_name = run_file_name_tokens{1}{3};
    
    % Relative paths used as unique identifiers for later naming of other intermediate and results files:
    run_relative_path = [run_file_name, run_file_extension];
    run_relative_path = clean_relative_path(run_relative_path);
    run_relative_path = regexprep(run_relative_path, '/+$', '');
    run_relative_path_list{run_index} = run_relative_path;
    run_key = strrep(run_relative_path, '/', '_');
    run_key_list{run_index} = run_key;

    % Some worksheets/run names are date, protein, run:
    % Accomodate spaces in the date:
    run_tokens = regexpi(run_file, '(([0-9]{2} ?){3}) *[^ ]* *run *([0-9]+) *\.csv', 'tokens');
    run_date = run_tokens{1}{1};
    run_date = strrep(run_date, ' ', '');
    run_number = run_tokens{1}{2};
    if master_script_options.debug_synapse_file_processing_verbose
      fprintf('Worksheet key "%s"\n', run_key)
    end

    if length(run_date) ~= 6
      warning('    run_date "%s" is not six digits, ignoring run', run_date)
      continue
    end
          
    % Check for run images in directories named after the protein or its synonyms:
    
    run_protein_name_synonyms_protein_index = find(cellfun(@(given_protein_name_synonyms)any(strcmp(given_protein_name_synonyms, run_protein)), master_script_options.protein_name_synonyms));
    if isempty(run_protein_name_synonyms_protein_index)
      run_protein_name_synonyms = {run_protein};
    else
      run_protein_name_synonyms = master_script_options.protein_name_synonyms{run_protein_name_synonyms_protein_index};
    end
    run_protein_name_number_synonyms = length(run_protein_name_synonyms);
    run_protein_name_synonym_index = 1;
    done = false;
    for image_location_index = 1:length(image_location)
      current_image_location = image_location{image_location_index};
      for run_protein_name_synonym_index = 1:run_protein_name_number_synonyms
        current_protein_synonym = run_protein_name_synonyms{run_protein_name_synonym_index};
        run_base_directory = clean_absolute_path([current_image_location, run_base_subdirectory, '/', current_protein_synonym, '/']);
        % done = exist(run_base_directory, 'file') > 0;
        if ~exist(run_base_directory, 'file')
          if master_script_options.debug_synapse_file_processing_verbose
            fprintf('    "%s" does not exist!\n', run_base_directory)
          end
          continue
        end
        
        [~, ~, protein_file_paths] = dirr(run_base_directory, 'name', '');
        protein_file_paths = protein_file_paths';
        protein_file_paths = strrep(protein_file_paths, run_base_directory, '');
        run_file_paths = protein_file_paths;
        run_paths_good = true(size(run_file_paths));
        % run_path_pattern = ['/[^/]*', run_date, '.*Run *', run_number, '[^/]*/'];
        run_path_pattern = sprintf('/[^/]*%s ?%s ?%s.*run *%s[^/]*/', run_date(1:2), run_date(3:4), run_date(5:6), run_number);
        run_paths_good = run_paths_good & ~cellfun(@isempty, regexpi(run_file_paths, run_path_pattern, 'start'));
        % run_file_paths
        run_file_paths = run_file_paths(run_paths_good);
        run_paths_good = true(size(run_file_paths));
        % run_file_paths
        
        % Find all fluorescence images of sensors:
        
        % image_path_patterns = {'/[^/]*T([0-9]{5})C[0-9]{2}Z\.[^/]*$', 'EGFP.*T([0-9]+)[^/]*\.[^/]+$', 'GFP.*\(Timepoint\W*([0-9]+)\)[^/]*\.[^/]+$'};
        image_path_patterns = {'/[^/]*T([0-9]{5})C[0-9]{2}Z\.[^/]*$', 'EGFP.*T([0-9]+)[^/]*\.[^/]+$', 'GFP.*\(Timepoint\W*([0-9]+)\)[^/]*\.[^/]+$'};

        for image_path_pattern_index = 1:length(image_path_patterns)
          image_path_pattern = image_path_patterns{image_path_pattern_index};
          run_image_paths_good = run_paths_good;
          run_image_paths_good = run_image_paths_good & ~cellfun(@isempty, regexp(run_file_paths, image_path_pattern, 'start'));
          run_image_paths = run_file_paths(run_image_paths_good);
          if isempty(run_image_paths)
            % Try another pattern:
            % warning('    pattern "%s" unmatched', image_path_pattern)
            continue
          else
            % Found images matching this pattern:
            run_image_paths = strcat(run_base_directory, run_image_paths);
          end
          % run_images_list{run_index} = run_image_paths;
          
          run_image_frame_indices = cellfun(@(x)str2num(x{1}{1}), regexp(run_image_paths, image_path_pattern, 'tokens'));
          
          if run_image_frame_indices(1) == 0
            if master_script_options.debug_synapse_file_processing_verbose
              fprintf('    run_key "%s": First frame is frame zero, adding one to the frame index for each image\n', run_key)
            end
            run_image_frame_indices = run_image_frame_indices + 1;
          end
          
          % Correct indices in run_images_list{run_index} so that indices corresponds to the frame number in the filename (now an empty image filename will appear where the image is missing):
          % Note: run_image_frame_indices_list is currently not used in other files, so everywhere it is assumed that all frames are present!
          run_image_frame_indices_list{run_index} = run_image_frame_indices;
          listed_run_image_paths = run_image_paths;
          run_image_paths = cell(max(run_image_frame_indices), 1);
          run_image_paths(run_image_frame_indices) = listed_run_image_paths;
          run_images_list{run_index} = run_image_paths;
          
        end
        
        if isempty(run_images_list{run_index})
          continue
        end

        % Find all brightfield images:
        
        % brightfield_image_path_pattern = [run_date, '.+Run *', run_number, '.*DIC.*/[^/]*\.[^/]*$'];
        % brightfield_image_path_pattern = [run_date, '.+Run *', run_number, '.*DIC.*/[^/]*DIC[^/]*\.[^/]*$'];
        % brightfield_image_path_pattern = [run_date, '.+Run *', run_number, '.*(DIC.*/)?[^/]*DIC\.[^/]*$'];
        brightfield_image_path_pattern = ['.*(DIC.*/)?[^/]*DIC\.[^/]*$'];
        run_brightfield_image_paths_good = run_paths_good;
        run_brightfield_image_paths_good = run_brightfield_image_paths_good & ~cellfun(@isempty, regexpi(run_file_paths, brightfield_image_path_pattern, 'start'));
        run_brightfield_image_paths = run_file_paths(run_brightfield_image_paths_good);
        if ~isempty(run_brightfield_image_paths)
          run_brightfield_image_paths = strcat(run_base_directory, run_brightfield_image_paths);
        end
        run_brightfield_images_list{run_index} = run_brightfield_image_paths;
        
        % Get the synapse locations in 2D and time:
        try
          [~, ~, run_file_extension] = fileparts(run_file);
          % [~, ~, run_data] = xlsread(run_file.file, run_file.worksheet, '', 'basic');
          % run_data = importdata(run_file, ',');
          % run_data = run_data.textdata;
          % run_data = cellfun(@conditional_num2str, run_data, 'UniformOutput', false);
          % run_data_number_columns_inference_data = importdata(run_file, ',');
          % if isfield(run_data_number_columns_inference_data, 'colheaders')
            % run_data_number_columns_inference_data = getfield(run_data_number_columns_inference_data, 'colheaders');
          % else
            % run_data_number_columns_inference_data = getfield(run_data_number_columns_inference_data, 'textdata');
          % end
          % run_data_number_columns = size(run_data_number_columns_inference_data, 2);
          run_file_handle = fopen(run_file);
          run_data_first_line = fgetl(run_file_handle);
          % Find number of columns as number of non-quoted commas plus one:
          run_data_first_line_commas = regexprep(run_data_first_line, '"[^"]*"', '');
          % run_data_first_line_commas = regexprep(run_data_first_line, '("[^"]*")|(''[^"]*'')', '');
          % run_data_first_line_commas = regexprep(run_data_first_line_commas, '[^,]+', '');
          % beep, keyboard
          % run_data_number_columns = length(run_data_first_line_commas) + 1;
          run_data_number_columns = sum(run_data_first_line_commas == ',') + 1;
          fseek(run_file_handle, 0, 'bof');
          % run_data = textscan(run_file_handle, [repmat('%q,', 1, run_data_number_columns - 1), '%s']);
          % run_data = textscan(run_file_handle, '%s', 'Delimiter', ',');
          % run_data = textscan(run_file_handle, [repmat('%[^,"\r\n]', 1, run_data_number_columns - 1), '%[^,"\r\n]'], inf);
          run_data = {};
          field_pattern = '([^,"\r\n]*|("[^"\r\n]*"))';
          while ~feof(run_file_handle)
            % run_data{end + 1} = textscan(run_file_handle, '%q', 'Delimiter', ',');
            current_line = fgetl(run_file_handle);
            % run_data{end + 1, 1} = textscan(current_line, [repmat('%[^,"\r\n],', 1, run_data_number_columns - 1), '%[^,"\r\n]']);
            run_data{end + 1, 1} = regexp(current_line, [repmat([field_pattern, ','], 1, run_data_number_columns - 1), field_pattern], 'tokens');
            run_data{end} = run_data{end}{1};
            % current_line
            % pause
            % beep, keyboard
          end
          fclose(run_file_handle);
          run_data = cat(1, run_data{:});
          % run_data = run_data.textdata;
          run_data_empty = cellfun(@isempty, run_data);
          % run_data(run_data_empty) = {{nan}};
          run_data(run_data_empty) = {nan};
          % run_data = cellfun(@(x)x{1}, run_data, 'UniformOutput', false);
          run_data = cellfun(@conditional_str2double, run_data, 'UniformOutput', false);
          % whos run_data, beep, keyboard
        catch err
          fprintf('    Error reading coordinates!\n')
          rethrow(err)
        end
        [run_data, run_data_format] = get_excel_file_coordinates(run_data);
        
        ignore_runs_without_brightfield_images = false;
        % ignore_runs_without_brightfield_images = true;
        if isempty(run_brightfield_images_list{run_index})
          if ignore_runs_without_brightfield_images
            warning('    No brightfield stack detected, ignoring run')
            continue
          else
            if master_script_options.debug_synapse_file_processing_verbose
              fprintf('    No brightfield stack detected\n')
            end
          end
        end

        
        done = true;
        
        % % Used to break once we found the directory, but now we might have multiple entries in image_location, so do not:
        if done, break, end
      end
      % % Used to break once we found the directory, but now we might have multiple entries in image_location, so do not:
      if done, break, end
    end

    if done
      runs_to_keep(run_index) = true;
    else
      if master_script_options.debug_synapse_file_processing_verbose
        warning('    No appropriate images found, ignoring run')
      end
    end
    
    % error('Implementation yet unfinished below this line!')
    
    
    run_data_list{run_index} = run_data;
    run_data_format_list{run_index} = run_data_format;
    run_data_as_array = cell2mat(reshape(run_data, [], 1));
    if isempty(run_data_as_array)
      if master_script_options.debug_synapse_file_processing_verbose
        % In case there are no cell couples, produce a warning:
        fprintf('    No cell couples specified, ignoring run\n')
      end
      runs_to_keep(run_index) = false;
      continue
    end
    
    % Filter cells with one frame:
    if master_script_options.minimum_number_frames_per_cell > 1
      number_frames_per_cell = hist(run_data_as_array(:, end), unique(run_data_as_array(:, end)));
      cells_to_remove = reshape(find(number_frames_per_cell < master_script_options.minimum_number_frames_per_cell), 1, []);
      for cell_to_remove = cells_to_remove
        for frame_index = 1:length(run_data)
          if numel(run_data{frame_index}) == 0
            continue
          end
          indices_to_remove = run_data{frame_index}(:, end) == cell_to_remove;
          run_data{frame_index}(indices_to_remove, :) = [];
        end
      end
      run_data_list{run_index} = run_data;
      run_data_as_array = cell2mat(reshape(run_data, [], 1));
    end
    
    if isempty(run_data_as_array)
      % In case there are no cell couples left, produce a warning:
      warning('    No cell couples specified after filtering cells with fewer than %d frames, ignoring run', master_script_options.minimum_number_frames_per_cell)
      runs_to_keep(run_index) = false;
      continue
    end
    
    number_cells = max(run_data_as_array(:, end));
    number_frames = size(run_data, 1);
    run_number_cells(run_index) = number_cells;
    run_number_frames(run_index) = number_frames;
    
    % Make sure the image list has a blank entry for every frame entry in run_data:
    if length(run_images_list{run_index}) < number_frames
      run_images_list{run_index} = [run_images_list{run_index}; cell(number_frames - size(run_images_list{run_index}, 1), 1)];
    end

    
    % This assumes the files are named with numbers 1 to number_frames:
    cell_synapse_point_tracks = cell(number_cells, number_frames);
    cell_relative_time_tracks = cell(number_cells, number_frames);
    
    % Produce a text file listing the position of each cell in its first frame:
    cell_indices_seen = [];
    % X, Y, frame number
    cell_positions = zeros(0, 3);
    
    
    % Loop through all frames with marked synapses and assign them to one time series per cell:
    for frame_index = 1:number_frames
      frame_synapse_locations = run_data{frame_index};
      for frame_cell_index = 1:size(frame_synapse_locations, 1)
      
        % For any method listing points, the last two are current_cell_relative_time and current_cell_index, the rest being the landmarks, with the first two assumed to be integer-valued 2D coordinates suitable for use with the old one-point annotation pipeline:
        current_cell_relative_time = frame_synapse_locations(frame_cell_index, end - 1);
        if isempty(find(timepoints_to_include == current_cell_relative_time))
            continue;
        end       
        current_cell_index = frame_synapse_locations(frame_cell_index, end);
        % Add segmentation for this frame to the appropriate cell:
        current_synapse_center_rounded = frame_synapse_locations(frame_cell_index, 1:2);
        cell_relative_time_tracks{current_cell_index, frame_index} = current_cell_relative_time;
        
        if master_script_options.use_two_point_synapses
          cell_synapse_point_tracks{current_cell_index, frame_index} = frame_synapse_locations(frame_cell_index, 1:end - 2);
        else
          cell_synapse_point_tracks{current_cell_index, frame_index} = current_synapse_center_rounded;
        end

        % Produce a text file listing the position of each cell in its first frame:
        if master_script_options.create_cell_position_text_file
          if ~ismember(current_cell_index, cell_indices_seen)
            cell_indices_seen = [cell_indices_seen; current_cell_index];
            cell_positions(end + 1, :) = [current_synapse_center_rounded, frame_index];
          end
        end
        
      end
    end

    if master_script_options.create_cell_position_text_file
      cell_positions_filename = [t_cell_info.path_info.cell_positions_location, run_key_list{run_index}, '.csv'];
      file_handle = fopen(cell_positions_filename, 'wt');
      fprintf(file_handle, 'X, Y, frame number\n');
      fclose(file_handle);
      dlmwrite(cell_positions_filename, cell_positions, '-append');
    end
    
    run_point_tracks{run_index} = cell_synapse_point_tracks;
    run_relative_time_tracks{run_index} = cell_relative_time_tracks;
    
    % How many cells are available for each relative time?
    run_relative_times{run_index} = unique(cell2mat(cell_relative_time_tracks(:)));
    if any(~ismember(run_relative_times{run_index}(:), [-2:9]))
      warning('run_relative_times{%d} (run_key "%s") contains non-standard relative times!', run_index, run_key_list{run_index})
    end
    run_number_relative_times(run_index) = length(run_relative_times{run_index});
    
    % What relative times do we have overall?
    all_relative_times = unique([all_relative_times; reshape(run_relative_times{run_index}, [], 1)]);
    number_all_relative_times = length(all_relative_times);
    
    % Collect relative times to determine approximately how many seconds each frame is from synapse formation:
    for relative_times_row_index = 1:size(cell_relative_time_tracks, 1)
      current_relative_times = cell2mat(cell_relative_time_tracks(relative_times_row_index, :));
      current_relative_times_zero_index = find(current_relative_times == 0);
      if isempty(current_relative_times_zero_index)
        continue
      end
      current_frames = find(~cellfun(@isempty, cell_relative_time_tracks(relative_times_row_index, :)));
      current_frames = (current_frames - current_frames(current_relative_times_zero_index)) * 20;

      all_relative_times_seconds_data = [all_relative_times_seconds_data; current_relative_times', current_frames'];
      
    end
    
  end
  
  all_relative_times_seconds = zeros(size(all_relative_times));
  for relative_time_index = 1:number_all_relative_times
    all_relative_times_seconds(relative_time_index) = mean(all_relative_times_seconds_data(all_relative_times_seconds_data(:, 1) == all_relative_times(relative_time_index), 2));
  end
  all_relative_times_seconds
  all_relative_times_seconds_rounded = round(all_relative_times_seconds)
  
  
  function filter_runs(given_runs_to_keep)
    % Remove runs from all relevant variables compactly:
    for filterable_run_variable_index = 1:length(filterable_run_variables)
      filtering_command = sprintf('%s = %s(given_runs_to_keep, :);', filterable_run_variables{filterable_run_variable_index * [1, 1]});
      eval(filtering_command);
    end
  end
  
  
  % Remove runs without images:
  
  filter_runs(runs_to_keep);
  
  number_runs = length(run_relative_path_list);
  
  
  % Filename and run_key parsing:
  
  condition_retrieval_function = @(given_strings)cellfun(@(x)x{1}{1}, regexp(given_strings, '[^_/]+ - ([^_/]+)_[^_/]+_[^_/]+\.csv', 'tokens'), 'UniformOutput', false);
  sensor_retrieval_function = @(given_strings)cellfun(@(x)x{1}{1}, regexp(given_strings, '[^_/]+ - [^_/]+_([^_/]+)_[^_/]+\.csv', 'tokens'), 'UniformOutput', false);
  % Assume everything was taken this century, and reorder so that datenum operates directly on this output after running cell2mat:
  % run_date_retrieval_function = @(given_strings)cellfun(@(x)cellfun(@(y)str2num(y), x{1}(1:3)) + [0, 0, 2000], regexp(given_strings, '[^_/]+ - [^_/]+_[^_/]+_[^_/]*([0-9]{2}) *([0-9]{2}) *([0-9]{2}) *[^_/]* *[Rr]un *[0-9]+ *\.csv', 'tokens'), 'UniformOutput', false);
  run_date_retrieval_function = @(given_strings)cellfun(@(x)cellfun(@(y)str2num(y), x{1}([3, 1, 2])) + [2000, 0, 0], regexp(given_strings, '[^_/]+ - [^_/]+_[^_/]+_[^_/]*([0-9]{2}) *([0-9]{2}) *([0-9]{2}) *[^_/]* *[Rr]un *[0-9]+ *\.csv', 'tokens'), 'UniformOutput', false);
  
  % Observed combinations, i.e.:
  % condition_sensor_combinations = unique(deformations_aligned_condition_sensor_list)
  % This doesn't work for some stupid reason:
  % condition_sensor_combinations = unique([condition_retrieval_function(run_key_list), sensor_retrieval_function(run_key_list)], 'rows')
  
  function [result_combinations] = unique_condition_sensor_combination_function(given_combinations)
    % Find unique condition-sensor combinations as represented by a two-column cell array of strings (first column condition, second sensor):
    % Explicit two-column version:
    result_combinations = strcat(given_combinations(:, 1), '|', given_combinations(:, 2));
    result_combinations = unique(result_combinations);
    [result_conditions, result_sensors] = strtok(result_combinations, '|');
    result_sensors = strrep(result_sensors, '|', '');
    result_combinations = [result_conditions, result_sensors];
  end
  
  function [result_combination_pairs] = unique_condition_sensor_combination_pair_function(given_combination_pairs)
    % % Find unique condition-sensor combinations as represented by a two-column cell array of strings (first column condition, second sensor):
    % Explicit two-by-two version:
    result_combination_pairs = cellfun(@(x)reshape(x.', 1, []), given_combination_pairs, 'UniformOutput', false);
    result_combination_pairs = cellfun(@(x)strjoin(x, '|'), result_combination_pairs, 'UniformOutput', false);
    result_combination_pairs = unique(result_combination_pairs);
    result_combination_pairs = cellfun(@(x)strsplit(x, '|'), result_combination_pairs, 'UniformOutput', false);
    result_combination_pairs = cellfun(@(x)reshape(x, 2, 2).', result_combination_pairs, 'UniformOutput', false);
  end
  
  function [result_memberships] = condition_sensor_combination_ismember(given_combinations1, given_combinations2)
    % Check if each combination from given_combinations1 is in given_combinations2 (assumes both are cell arrays with two columns):
    % % If the first entry is just a pair, assume everything else is:
    conditions1 = given_combinations1(:, 1);
    sensors1 = given_combinations1(:, 2);
    conditions2 = given_combinations2(:, 1);
    sensors2 = given_combinations2(:, 2);
    result_memberships = ismember(conditions1, conditions2) & ismember(sensors1, sensors2);
  end
  
  function [result_memberships] = condition_sensor_combination_pair_ismember(given_combination_pairs1, given_combination_pairs2)
    % Check if each pair of combinations from given_combination_pairs1 is in given_combination_pairs2 (assumes both are cell arrays of cell arrays):
    % Explicit two-by-two version:
    % % If the first entry is just a pair, assume everything else is:
    former_conditions1 = cellfun(@(x)x{1, 1}, given_combination_pairs1, 'UniformOutput', false);
    former_sensors1 = cellfun(@(x)x{1, 2}, given_combination_pairs1, 'UniformOutput', false);
    latter_conditions1 = cellfun(@(x)x{2, 1}, given_combination_pairs1, 'UniformOutput', false);
    latter_sensors1 = cellfun(@(x)x{2, 2}, given_combination_pairs1, 'UniformOutput', false);
    former_conditions2 = cellfun(@(x)x{1, 1}, given_combination_pairs2, 'UniformOutput', false);
    former_sensors2 = cellfun(@(x)x{1, 2}, given_combination_pairs2, 'UniformOutput', false);
    latter_conditions2 = cellfun(@(x)x{2, 1}, given_combination_pairs2, 'UniformOutput', false);
    latter_sensors2 = cellfun(@(x)x{2, 2}, given_combination_pairs2, 'UniformOutput', false);
    result_memberships = ismember(former_conditions1, former_conditions2) & ismember(former_sensors1, former_sensors2) & ismember(latter_conditions1, latter_conditions2) & ismember(latter_sensors1, latter_sensors2);
  end
 
  condition_sensor_combinations = unique_condition_sensor_combination_function([condition_retrieval_function(run_key_list), sensor_retrieval_function(run_key_list)]);
  % error
  
  % Remove runs for specific sensors:
  
  run_keys_to_exclude = master_script_options.run_keys_to_exclude;
  for sensors_to_exclude_index = 1:length(sensors_to_exclude)
    fprintf('Removing sensor %s!\n', sensors_to_exclude{sensors_to_exclude_index})
    run_keys_to_exclude = [run_keys_to_exclude; run_key_list(strcmp(sensor_retrieval_function(run_key_list), sensors_to_exclude{sensors_to_exclude_index}))];
    condition_sensor_combinations(strcmp(condition_sensor_combinations(:, 2), sensors_to_exclude{sensors_to_exclude_index}), :) = [];
  end
  
  % Remove specific runs:
  
  runs_to_keep = true(number_runs, 1);
  for run_keys_to_exclude_index = 1:length(run_keys_to_exclude)
    fprintf('Removing run_key %s!\n', run_keys_to_exclude{run_keys_to_exclude_index})
    runs_to_keep(strcmp(run_key_list, run_keys_to_exclude{run_keys_to_exclude_index})) = false;
  end
  
  % Remove runs in specific date ranges:
  
  date_ranges_to_exclude = master_script_options.date_ranges_to_exclude;
  current_run_dates = run_date_retrieval_function(run_key_list);
  
  if any(cellfun(@isempty, current_run_dates))
    current_run_dates
    error('current_run_dates empty for at least one run!')
  end
  current_run_dates = cell2mat(current_run_dates);
  current_run_datenums = datenum(current_run_dates);
  date_range_runs_to_keep = true(number_runs, 1);
  for date_ranges_to_exclude_index = 1:size(date_ranges_to_exclude, 1)
    fprintf('Removing date_range %04d-%02d-%02d to %04d-%02d-%02d!\n', date_ranges_to_exclude(date_ranges_to_exclude_index, :))
    % runs_to_keep((current_run_datenums >= datenum(date_ranges_to_exclude(date_ranges_to_exclude_index, 1:3))) & (current_run_datenums <= datenum(date_ranges_to_exclude(date_ranges_to_exclude_index, 4:6)))) = false;
    % runs_to_keep((current_run_datenums >= datenum(date_ranges_to_exclude(date_ranges_to_exclude_index, 1:3))) & (current_run_datenums <= datenum(date_ranges_to_exclude(date_ranges_to_exclude_index, 4:6)))) = false;
    date_range_runs_to_keep((current_run_datenums >= datenum(date_ranges_to_exclude(date_ranges_to_exclude_index, 1:3))) & (current_run_datenums <= datenum(date_ranges_to_exclude(date_ranges_to_exclude_index, 4:6)))) = false;
  end
  runs_to_keep = runs_to_keep & date_range_runs_to_keep;

  number_condition_sensors = size(condition_sensor_combinations, 1);
  condition_sensor_combinations, number_condition_sensors
  
  % Remove runs from all relevant variables compactly:
  number_runs, number_runs_to_keep = sum(runs_to_keep)
  filter_runs(runs_to_keep);
  number_runs = number_runs_to_keep;
  
  
  conditions = unique(condition_sensor_combinations(:, 1));
  sensors = unique(condition_sensor_combinations(:, 2));
  number_conditions = length(conditions)
  number_sensors = length(sensors)
  
  % % Correct indices in run_images_list{run_index} so that indices corresponds to the frame number in the filename:
  
  fprintf('\n')
  
  % Get fluorescent image frame and voxel sizes:
  
  % warning('>>>> HACK, hard-coding images'' voxel sizes using the date of CW''s move to University of Bristol!')
  for run_index = 1:number_runs
  
    current_run_key = run_key_list{run_index};
    current_run_number_frames = run_number_frames(run_index);
    current_run_images_list = run_images_list{run_index};
    current_run_images_list = current_run_images_list(1:current_run_number_frames);
    current_run_brightfield_images_list = run_brightfield_images_list{run_index};
    current_run_brightfield_images_list = current_run_brightfield_images_list(1:min(end, current_run_number_frames));
    
    frames_available = ~cellfun(@isempty, current_run_images_list);
    fluorescent_image_frame = robust_read_stack(current_run_images_list{find(frames_available, 1)});
    fluorescent_image_frame_info = imfinfo(current_run_images_list{find(frames_available, 1)});
    % xruan 07/26/2015
    % disable the image frame description because in some downsampled
    % images, there is no such description ImageDescription. Maybe we can
    % add it later if possible. 
    % fluorescent_image_frame_description = strrep(fluorescent_image_frame_info(1).ImageDescription, char(13), '  ');
    fluorescent_image_frame_size = size(fluorescent_image_frame);
    
    run_image_size_list(run_index, :) = fluorescent_image_frame_size;
    current_run_date = cell2mat(run_date_retrieval_function({current_run_key}));
    % current_run_date(current_run_date(:, 3) < 100, 3) = current_run_date(current_run_date(:, 3) < 100, 3) + 2000;
    % All worksheet names' dates should be month day year (American style) as with old data:
    % current_run_date_number = datenum(current_run_date(3), current_run_date(1), current_run_date(2));
    current_run_date_number = datenum(current_run_date);
    % fprintf('current_run_date_number: %04d-%02d-%02d %d <= 2012-08-15 %d?\n', current_run_date([3, 1, 2]), current_run_date_number, datenum(2012, 8, 15))
    if current_run_date_number < datenum(2012, 8, 15) || current_run_date_number >= datenum(2014, 3, 9)
      run_image_voxel_size_list(run_index, :) = [0.406, 0.406, 0.4];
    else
      run_image_voxel_size_list(run_index, :) = [0.34, 0.34, 1];
    end
    % warning('2014-11-19: Add rule for HS1 images! Or delete associated .csv files!')
    
  end
  
  % Segmentation and such were designed for this size, scale everything accordingly:
  design_voxel_size = [0.406, 0.406, 0.4];
  % % We want to rescale images to have maximum resolution and cubical voxels, so take the minimum voxel dimension:
  % desired_voxel_size = min(run_image_voxel_size_list(:)) .* ones(1, 3);
  desired_voxel_size = design_voxel_size;  
  
  % Filter synapse centers outside of the image frame:
  
  runs_to_keep = true(number_runs, 1);
  
  for run_index = 1:number_runs
    
    current_run_key = run_key_list{run_index};
    current_run_number_frames = run_number_frames(run_index);
    current_run_images_list = run_images_list{run_index};
    current_run_images_list = current_run_images_list(1:current_run_number_frames);
    current_run_brightfield_images_list = run_brightfield_images_list{run_index};
    current_run_brightfield_images_list = current_run_brightfield_images_list(1:min(end, current_run_number_frames));
    current_run_point_tracks = run_point_tracks{run_index};
    
    frames_available = ~cellfun(@isempty, current_run_images_list);

    frame_size = run_image_size_list(run_index, :);
    
    frames_annotated = any(~cellfun(@isempty, current_run_point_tracks), 1).';
    cells_annotated = any(~cellfun(@isempty, current_run_point_tracks), 2);
    frame_cells_within_bounds = cellfun(@(x)(~isempty(x)) && (x(1, 1) >= 1 && x(1, 2) >= 1 && x(1, 1) <= frame_size(2) && x(1, 2) <= frame_size(1)), current_run_point_tracks);    current_run_point_tracks(~frame_cells_within_bounds) = {[]};
    run_point_tracks{run_index} = current_run_point_tracks;  
  end
  
  filter_runs(runs_to_keep);
  number_runs = number_runs_to_keep;
  
  % 02/25/2016 xruan make all used images as a single cell array and without including the
  % the empty cells.
  run_key_index_list = [];
  run_key_use_list = {};
  image_name_run_list = {};
  synapse_tracks_run_list = {};
  relative_time_run_list = {};
  image_frame_run_list = [];
  cell_number_run_list = [];
  
  for run_index = 1:number_runs  
    current_run_key = run_key_list{run_index};
    current_run_relative_time_tracks = run_relative_time_tracks{run_index};
    current_run_point_tracks = run_point_tracks{run_index};
    for time_ind = 1 : numel(timepoints_to_include)
        curr_time = timepoints_to_include(time_ind);
        frame_inds = cellfun(@(x) ~isempty(x) && x(1) == curr_time, current_run_relative_time_tracks);
        frame_inds_num = sum(frame_inds(:));
        [cell_num, frame_num] = ind2sub(size(frame_inds), find(frame_inds == 1));
        
        run_key_use_list(end + 1 : end + frame_inds_num, 1) = {current_run_key};
        image_name_run_list(end + 1 : end + frame_inds_num, 1) = run_images_list{run_index}(frame_num);
        synapse_tracks_run_list(end + 1 : end + frame_inds_num, 1) = run_point_tracks{run_index}(frame_inds == 1);
        relative_time_run_list(end + 1 : end + frame_inds_num, 1) = run_relative_time_tracks{run_index}(frame_inds == 1);
        image_frame_run_list  = [image_frame_run_list; frame_num];
        cell_number_run_list = [cell_number_run_list; cell_num];
        run_key_index_list = [run_key_index_list; run_index * ones(numel(cell_num), 1)];
    end 
  end
  
  
  t_cell_info.synapse_info.run_relative_path_list = run_relative_path_list;
  t_cell_info.synapse_info.run_key_list = run_key_list;
  t_cell_info.synapse_info.run_data_list = run_data_list;
  t_cell_info.synapse_info.run_images_list = run_images_list;
  t_cell_info.synapse_info.run_brightfield_images_list = run_brightfield_images_list;
  t_cell_info.synapse_info.run_image_frame_indices_list = run_image_frame_indices_list;
  t_cell_info.synapse_info.run_image_size_list = run_image_size_list;
  t_cell_info.synapse_info.run_image_voxel_size_list = run_image_voxel_size_list;
  t_cell_info.synapse_info.desired_voxel_size = desired_voxel_size;
  t_cell_info.synapse_info.design_voxel_size = design_voxel_size;
  
  t_cell_info.synapse_info.run_number_cells = run_number_cells;
  t_cell_info.synapse_info.run_number_frames = run_number_frames;
  t_cell_info.synapse_info.run_point_tracks = run_point_tracks;
  t_cell_info.synapse_info.run_relative_time_tracks = run_relative_time_tracks;
  t_cell_info.synapse_info.run_relative_times = run_relative_times;
  t_cell_info.synapse_info.run_number_relative_times = run_number_relative_times;
  t_cell_info.synapse_info.all_relative_times = all_relative_times;
  t_cell_info.synapse_info.number_all_relative_times = number_all_relative_times;
  t_cell_info.synapse_info.all_relative_times_seconds_data = all_relative_times_seconds_data;
  t_cell_info.synapse_info.all_relative_times_seconds = all_relative_times_seconds;
  t_cell_info.synapse_info.all_relative_times_seconds_rounded = all_relative_times_seconds_rounded;
  t_cell_info.synapse_info.number_runs = number_runs;
  
  t_cell_info.synapse_info.condition_sensor_combinations = condition_sensor_combinations;
  t_cell_info.synapse_info.number_condition_sensors = number_condition_sensors;
  t_cell_info.synapse_info.conditions = conditions;
  t_cell_info.synapse_info.sensors = sensors;
  t_cell_info.synapse_info.number_conditions = number_conditions;
  t_cell_info.synapse_info.number_sensors = number_sensors;
  t_cell_info.synapse_info.condition_retrieval_function = condition_retrieval_function;
  t_cell_info.synapse_info.sensor_retrieval_function = sensor_retrieval_function;
  t_cell_info.synapse_info.run_date_retrieval_function = run_date_retrieval_function;
  t_cell_info.synapse_info.unique_condition_sensor_combination_function = @unique_condition_sensor_combination_function;
  t_cell_info.synapse_info.unique_condition_sensor_combination_pair_function = @unique_condition_sensor_combination_pair_function;
  t_cell_info.synapse_info.condition_sensor_combination_ismember = @condition_sensor_combination_ismember;
  t_cell_info.synapse_info.condition_sensor_combination_pair_ismember = @condition_sensor_combination_pair_ismember;

  t_cell_info.synapse_info.run_key_use_list = run_key_use_list;
  t_cell_info.synapse_info.image_name_run_list = image_name_run_list;
  t_cell_info.synapse_info.synapse_tracks_run_list = synapse_tracks_run_list;
  t_cell_info.synapse_info.relative_time_run_list = relative_time_run_list;
  t_cell_info.synapse_info.image_frame_run_list = image_frame_run_list;
  t_cell_info.synapse_info.cell_number_run_list = cell_number_run_list;
  t_cell_info.synapse_info.run_key_index_list = run_key_index_list;
  
end  

