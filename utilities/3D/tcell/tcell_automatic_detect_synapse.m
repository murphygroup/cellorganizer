function [t_cell_info] = tcell_automatic_detect_synapse(t_cell_info, options)
  %  automatically detect synapse position, given one
  % synapse and relative time currently. 
  % 2018-08-20 xruan: copied from 
  % /projects/cellorganizer/xruan/tcell_project/tcell_clean_dev/master_script_automatic_detect_synapse.m
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
  inferred_synapse_location = t_cell_info.path_info.inferred_synapse_location;
  temp_location = t_cell_info.path_info.temp_location;
  
  sensors_to_exclude = master_script_options.sensors_to_exclude;
  
  t_cell_info.synapse_info = struct();
    
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
  
  run_folder_list = cell(number_runs, 1);

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
  chunk_lock_clean(inferred_synapse_location, 1 * 60);
  run_file_uncomputed = true(number_runs, 1);
  
  while any(run_file_uncomputed)
  
      for run_index = 1:number_runs
        % Filename of file specifying synapse locations in 2D and time:
        run_file = run_list{run_index};
        [~, run_filename] = fileparts(run_file);
        [can_start, final_name, final_exists] = chunk_start_clean(inferred_synapse_location, run_filename, '.csv');
        
        if final_exists
            run_file_uncomputed(run_index) = false;
            continue
        end

        if ~can_start
          continue
        end
        
        fprintf('Infer synapse info: %s\n', run_file)
        
        [~, run_file_name, run_file_extension] = fileparts(run_file);

        run_relative_path = [run_file_name, run_file_extension];
        run_relative_path = clean_relative_path(run_relative_path);
        run_relative_path = regexprep(run_relative_path, '/+$', '');
        run_relative_path_list{run_index} = run_relative_path;
        run_key = strrep(run_relative_path, '/', '_');
        run_key_list{run_index} = run_key;
    
        cur_run_result_location = [inferred_synapse_location, '/', run_file_name, '/'];
        mkdir(cur_run_result_location);
        % Some worksheets/run names are date, protein, run:
        % Accomodate spaces in the date:
        
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
        
        run_data_options = struct();
        run_data_options.run_index = run_index;
        run_data_options.run_file = run_file;
        [run_data, run_data_by_category, run_data_format] = get_excel_file_coordinates_new(run_data, run_data_options);
        
        % find the run data info for starting time points, and convert the
        % input filename into the format for the pipeline. 
        starting_time_point = 0;
        start_points_filename = [cur_run_result_location, 'start_', run_file_name, '.csv'];
        [start_points_filename, run_folder] = convert_run_data_to_pipeline_csv_input_format(run_data, starting_time_point, start_points_filename); 
                
        done = false;
        run_folder_list{run_index} = run_folder;

        [~, filename] = fileparts(run_file);
        annotation_shortname = [filename, '.csv'];
        if master_script_options.infer_starting_time_point
            start_points_filename = [temp_location, 'start_', annotation_shortname];
            synapse_result_filename = [inferred_synapse_location, annotation_shortname];
            if ~exist(start_points_filename, 'file')
                detect_starting_point(run_folder, start_points_filename)
            end 
            detect_param.verbose = false;
            detect_param.intermediate_location = temp_location;
            detect_param.save_intermediate_results = false;
            detect_param.input_mode = 'new';
            detect_param.filter_touching_tcell = true;
            [inferred_csv_filename] = infer_synapse_location(start_points_filename, run_folder, temp_location, detect_param);
        else             
            % [annotation_name] = convert_csv_format_to_pipeline(run_file, temp_location);
            % start_points_filename = [temp_location, annotation_name];
            synapse_result_filename = [inferred_synapse_location, annotation_shortname];
            detect_param.verbose = false;
            detect_param.intermediate_location = cur_run_result_location;
            detect_param.save_intermediate_results = false;
            detect_param.input_mode = 'new';
            detect_param.filter_touching_tcell = true;
            [inferred_csv_filename] = infer_synapse_location(start_points_filename, run_folder, cur_run_result_location, detect_param);
        end
        [reformat_csv_filename] = convert_inferred_synapse_file_to_cellorganizer_format(inferred_csv_filename, run_folder, synapse_result_filename);
        % [errors] = calculate_inferred_synapse_difference_3(run_file, synapse_result_filename);               
            
        if false
            predict_csv_filename = synapse_result_filename;
            raw_annotation_filename = run_file;
            run_folder = run_folder;
            temp_figure_dir = [temp_location, 'figures/'];
            mkdir(temp_figure_dir);
            temp_figure_original_dir = [temp_location, 'figures_original/'];
            mkdir(temp_figure_original_dir);
            [~, filename] = fileparts(raw_annotation_filename);
            save_folder = [temp_figure_dir, filename, '/'];
            mkdir(save_folder);
            save_original_folder = [temp_figure_original_dir, filename, '/'];
            mkdir(save_original_folder);
            if master_script_options.infer_starting_time_point
                % plot_cell_shape_change_measurement(predict_csv_filename, raw_annotation_filename, run_folder, save_folder);                     
                visualize_synapse_prediction_1(predict_csv_filename, run_folder, save_folder);                     
                visualize_synapse_prediction_1(raw_annotation_filename, run_folder, save_original_folder);                     
                visualize_synapse_prediction_3(start_points_filename, predict_csv_filename, raw_annotation_filename, run_folder, save_original_folder);                     
            else
                visualize_synapse_prediction(predict_csv_filename, raw_annotation_filename, run_folder, save_folder);                    
            end
        end

        done = true;        

        if false && strcmpi(run_protein, 'Actin')
          % Debug info:
          % run_base_directory
          run_file.file, run_file.worksheet, run_base_directory, done
          keyboard
        end

      if done
          runs_to_keep(run_index) = true;
          run_file_uncomputed(run_index) = false;
      else
          if master_script_options.debug_synapse_file_processing_verbose
            warning('    No appropriate images found, ignoring run')
          end
      end

      % error('Implementation yet unfinished below this line!')
      % save(final_name, 'errors', 'run_file');
      chunk_finish(inferred_synapse_location, run_filename);

    end
    if any(run_file_uncomputed)
      out_status = chunk_lock_clean(inferred_synapse_location, 20);
      disp('Wait for other running to fininsh!');
      pause(20);
    end
  end
  % error('Below has not been implemented!');
  
  if false
    % Debug info:
    % run_data
    run_data_list
    run_data_format_list
    run_point_tracks
    run_relative_time_tracks
    beep, keyboard
  end
  
  
end  


function [start_points_filename, run_folder] = convert_run_data_to_pipeline_csv_input_format(run_data, starting_time_point, start_points_filename)
% convert the format of run data to just save information for start time
% points, and also extract run folder. 
% by default we assume in a folder there is only one movie, that is, all
% movie should be within the same directory for a single run data file. 

run_data_relative_time_mat = [run_data{:, 4}];

run_data_start_time = run_data(run_data_relative_time_mat == starting_time_point, :);

% extract run folder
run_data_filenames = run_data_start_time(:, 1);
run_folder_cell = cellfun(@fileparts, run_data_filenames, 'UniformOutput', false);

run_folder = unique(run_folder_cell);
if numel(run_folder) > 1 
    error('A run file can only contain one movie!');
end

% remove the GFP level. 
run_folder = [fileparts(run_folder{1}), '/'];

% find frame numbers for all starting time points by regex matching. 
frame_names = regexp(run_data_filenames, '([^/]*$)', 'tokens');
frame_names = cellfun(@(x) x{1}{1}, frame_names, 'uniformoutput', false);

regexp_pattern = '([0-9]{1,3})';

frame_nums = regexp(frame_names, regexp_pattern, 'tokens');
frame_nums = cellfun(@(x) str2num(x{1}{1}), frame_nums);


% left and right end-points
end_point_coords = run_data_start_time(:, 3);
end_point_coords = cat(1, end_point_coords{:});
csv_content = [end_point_coords, frame_nums];

csv_title={'Lx','Ly','Rx','Ry','start_frame'};
csv_content=num2cell(csv_content);
output = cell2table(csv_content,'VariableNames',csv_title);
writetable(output, start_points_filename,'WriteRowNames',true);


end

