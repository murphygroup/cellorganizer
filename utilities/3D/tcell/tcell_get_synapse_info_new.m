function [t_cell_info] = tcell_get_synapse_info_new(t_cell_info, options)
  % 2016-02-25 xruan: Copied from master_script_get_synapse_info.m.
  % 2/2/2021 R.F. Murphy - allow timepoints_to_include to be any ("*")
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
  
  % 08/20/2018 set whether to use inferred synapse
  if master_script_options.infer_synapses
      original_synapse_location = synapse_location;
      inferred_synapse_location = [t_cell_info.path_info.inferred_synapse_location, '*.csv'];
      synapse_location = {inferred_synapse_location};
      t_cell_info.options.original_synapse_location = original_synapse_location;
      t_cell_info.options.synapse_location = synapse_location;
  end  
  
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
  
  t_cell_info.synapse_info.clean_relative_path = @clean_relative_path;
  t_cell_info.synapse_info.clean_absolute_path = @clean_absolute_path;  
  
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
  
  % Segmentation and such were designed for this size, scale everything accordingly:
  design_voxel_size = [0.406, 0.406, 0.4];
  % % We want to rescale images to have maximum resolution and cubical voxels, so take the minimum voxel dimension:
  % desired_voxel_size = min(run_image_voxel_size_list(:)) .* ones(1, 3);
  desired_voxel_size = design_voxel_size;  
  
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

    run_relative_path = [run_file_name, run_file_extension];
    run_relative_path = clean_relative_path(run_relative_path);
    run_relative_path = regexprep(run_relative_path, '/+$', '');
    run_relative_path_list{run_index} = run_relative_path;
    run_key = strrep(run_relative_path, '/', '_');
    run_key_list{run_index} = run_key;

    % Some worksheets/run names are date, protein, run:
    % Accomodate spaces in the date:

    if master_script_options.debug_synapse_file_processing_verbose
      fprintf('Worksheet key "%s"\n', run_key)
    end
          
    % Check for run images in directories named after the protein or its synonyms:
    
    done = false;
        
    % Get the synapse locations in 2D and time:
    try
      % [~, ~, run_file_extension] = fileparts(run_file);
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
    % [run_data, run_data_format] = get_excel_file_coordinates(run_data);
    run_data_options = struct();
    run_data_options.run_index = run_index;
    run_data_options.run_file = run_file;
    [run_data, run_data_by_category, run_data_format] = get_excel_file_coordinates_new(run_data, run_data_options);

    done = true;
        
      
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
    run_data_cell_index = cell2mat(run_data(:, 5));
    run_image_voxel_size_list(run_index, :) = desired_voxel_size;
    if isempty(run_data_cell_index)
      if master_script_options.debug_synapse_file_processing_verbose
        % In case there are no cell couples, produce a warning:
        fprintf('    No cell couples specified, ignoring run\n')
      end
      runs_to_keep(run_index) = false;
      continue
    end
    
    % Filter cells with one frame:
    % master_script_options.minimum_number_frames_per_cell = 12
    if master_script_options.minimum_number_frames_per_cell > 1
      number_frames_per_cell = hist(run_data_cell_index(:, end), unique(run_data_cell_index(:, end)));
      cells_to_remove = reshape(find(number_frames_per_cell < master_script_options.minimum_number_frames_per_cell), 1, []);
      frame_to_remove_indices = arrayfun(@(x) find(run_data_cell_index(:, end) == x), cells_to_remove, 'UniformOutput', false);
      frame_to_remove_indices = cell2mat(frame_to_remove_indices');
      run_data(frame_to_remove_indices, :) = [];
      run_data_cell_index(frame_to_remove_indices, :) = [];
    end
    
    if isempty(run_data_cell_index)
      % In case there are no cell couples left, produce a warning:
      warning('    No cell couples specified after filtering cells with fewer than %d frames, ignoring run', master_script_options.minimum_number_frames_per_cell)
      runs_to_keep(run_index) = false;
      continue
    end
    
    number_cells = numel(unique(run_data_cell_index(:, end)));
    number_frames = size(run_data, 1);
    run_number_cells(run_index) = number_cells;
    run_number_frames(run_index) = number_frames;
    
    % Make sure the image list has a blank entry for every frame entry in run_data:
    if length(run_images_list{run_index}) < number_frames
      run_images_list{run_index} = [run_images_list{run_index}; cell(number_frames - size(run_images_list{run_index}, 1), 1)];
    end    
    
    run_data_list{run_index} = run_data;

    cell_coordinates = cell2mat(run_data(:, 3));
    cell_synapse_center = [mean(cell_coordinates(:, [1, 3]), 2), mean(cell_coordinates(:, [2, 4]), 2)];
    cell_relative_time = cell2mat(run_data(:, 4));
    % X, Y, frame number
    cell_positions = [cell_synapse_center, cell_relative_time];
    if master_script_options.create_cell_position_text_file
      cell_positions_filename = [t_cell_info.path_info.cell_positions_location, run_key_list{run_index}, '.csv'];
      file_handle = fopen(cell_positions_filename, 'wt');
      fprintf(file_handle, 'X, Y, frame number\n');
      fclose(file_handle);
      dlmwrite(cell_positions_filename, cell_positions, '-append');
    end
  end
 
  
  % 02/25/2016 xruan make all used images as a single cell array and without including the
  % the empty cells.
  image_name_run_list = {};
  synapse_tracks_run_list = {};
  relative_time_run_list = {};
  image_frame_run_list = [];
  cell_number_run_list = [];
  
  % check if the setted timepoints_to_include exists, if not throw out
  % warning or error
  all_unqiue_timepoints = [];
  for run_index = 1:number_runs  
    run_data = run_data_list{run_index};
    frame_reltive_times = unique(cell2mat(run_data(:, 4)));
    all_unqiue_timepoints = unique([all_unqiue_timepoints; frame_reltive_times]);
  end
  if strcmp(timepoints_to_include,"any")
      timepoints_to_include = all_unqiue_timepoints;
  else
      non_exist_timepoints = setdiff(timepoints_to_include, all_unqiue_timepoints);
      if ~isempty(non_exist_timepoints)
          error('Time point %s not exist! Please check your code and set right time point(s)', mat2str(non_exist_timepoints));
      end
  end
  
  
  all_run_data = [];
  for run_index = 1:number_runs  
    current_run_key = run_key_list{run_index};
    run_data = run_data_list{run_index};
    frame_reltive_times = cell2mat(run_data(:, 4));
%    frame_inds_to_keep = sort(cell2mat(arrayfun(@(x) find(frame_reltive_times == x),  timepoints_to_include, 'UniformOutput', false)'));
    frame_inds_to_keep = [];
    for tti=1:length(timepoints_to_include)
        frame_inds_to_keep = [frame_inds_to_keep; ...
            find(frame_reltive_times==timepoints_to_include(tti))];
    end
    frame_inds_to_keep = sort(frame_inds_to_keep);
    frame_inds_num = numel(frame_inds_to_keep);
    cur_run_data = run_data(frame_inds_to_keep, :);
    all_run_data = [all_run_data; cur_run_data];
    
    run_key_list(end + 1 : end + frame_inds_num, 1) = {current_run_key};
%     image_name_run_list(end + 1 : end + frame_inds_num, 1) = cur_run_data(:, 1);
%     relative_time_run_list(end + 1 : end + frame_inds_num, 1) = run_relative_time_tracks{run_index}(frame_inds == 1);
%     image_frame_run_list  = [image_frame_run_list; frame_num];
%     cell_number_run_list = [cell_number_run_list; cell_num];
%     run_key_index_list = [run_key_index_list; run_index * ones(numel(cell_num), 1)];
  end
  image_name_run_list = all_run_data(:, 1);
  frame_channel_list = cell2mat(all_run_data(:, 2));
  synapse_tracks_run_list = cell2mat(all_run_data(:, 3));
  relative_time_run_list = cell2mat(all_run_data(:, 4));
  frame_tracking_list = cell2mat(all_run_data(:, 5));
  
  
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
  
%   t_cell_info.synapse_info.run_number_cells = run_number_cells;
%   t_cell_info.synapse_info.run_number_frames = run_number_frames;
%   t_cell_info.synapse_info.run_point_tracks = run_point_tracks;
%   t_cell_info.synapse_info.run_relative_time_tracks = run_relative_time_tracks;
%   t_cell_info.synapse_info.run_relative_times = run_relative_times;
%   t_cell_info.synapse_info.run_number_relative_times = run_number_relative_times;
%   t_cell_info.synapse_info.all_relative_times = all_relative_times;
%   t_cell_info.synapse_info.number_all_relative_times = number_all_relative_times;
%   t_cell_info.synapse_info.all_relative_times_seconds_data = all_relative_times_seconds_data;
%   t_cell_info.synapse_info.all_relative_times_seconds = all_relative_times_seconds;
%   t_cell_info.synapse_info.all_relative_times_seconds_rounded = all_relative_times_seconds_rounded;
%   t_cell_info.synapse_info.number_runs = number_runs;
  
%   t_cell_info.synapse_info.condition_sensor_combinations = condition_sensor_combinations;
%   t_cell_info.synapse_info.number_condition_sensors = number_condition_sensors;
%   t_cell_info.synapse_info.conditions = conditions;
%   t_cell_info.synapse_info.sensors = sensors;
%   t_cell_info.synapse_info.number_conditions = number_conditions;
%   t_cell_info.synapse_info.number_sensors = number_sensors;
%   t_cell_info.synapse_info.condition_retrieval_function = condition_retrieval_function;
%   t_cell_info.synapse_info.sensor_retrieval_function = sensor_retrieval_function;
%   t_cell_info.synapse_info.run_date_retrieval_function = run_date_retrieval_function;
%   t_cell_info.synapse_info.unique_condition_sensor_combination_function = @unique_condition_sensor_combination_function;
%   t_cell_info.synapse_info.unique_condition_sensor_combination_pair_function = @unique_condition_sensor_combination_pair_function;
%   t_cell_info.synapse_info.condition_sensor_combination_ismember = @condition_sensor_combination_ismember;
%   t_cell_info.synapse_info.condition_sensor_combination_pair_ismember = @condition_sensor_combination_pair_ismember;

  t_cell_info.synapse_info.all_run_data = all_run_data;
  t_cell_info.synapse_info.image_name_run_list = image_name_run_list;
  t_cell_info.synapse_info.frame_channel_list = frame_channel_list;  
  t_cell_info.synapse_info.synapse_tracks_run_list = synapse_tracks_run_list;
  t_cell_info.synapse_info.relative_time_run_list = relative_time_run_list;
  t_cell_info.synapse_info.frame_tracking_list = frame_tracking_list;
  % t_cell_info.synapse_info.run_key_index_list = run_key_index_list;
  t_cell_info.synapse_info.included_timepoints = timepoints_to_include;
  
end  


  