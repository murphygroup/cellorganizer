function [values, file_ids, finished_values] = traverse_directory_structure(options)
  % Apply the function options.compute_function to unmorphed images from the directory segmentations_path passing if should_compute_function on that image returns true. should_compute_function(file_id) takes one argument, a struct with fields date, protein, run, cell_number, frame_number, and filename. compute_function(segmentation_struct, file_id) takes the segmentation file's contents as the first argument and an identical file_id for an image passing the should_compute_function test.
  % 
  % File Exchange dependencies:
  % http://www.mathworks.com/matlabcentral/fileexchange/8682-dirr-find-files-recursively-filtering-name-date-or-bytes
  % http://www.mathworks.com/matlabcentral/fileexchange/31272-datahash
  % 
  % tebuck dependencies:
  % chunk_start.m, chunk_finish.m
  % 
  % To do:
  % Add check to see if loaded files are newer, if so, recompute, or unnecessary now?
  % 
  % 2013-01-29 tebuck: Copied from traverse_t_cell_segmented_images.m.
  % 2013-03-25 tebuck: Added options.save_filename_function.
  
  warning('off', 'MATLAB:MKDIR:DirectoryExists');

  % Process options:
  default_options = struct(); 
  default_options.load_location = 'temp1';
  default_options.save_location = 'temp2';
  default_options.load_pattern = 'temp2';
  default_options.just_save_results = false;
  default_options.just_run_finished_compute_function = false;
  default_options.dry_run = false;
  default_options.load_only = false;
  default_options.verbose = false;
  default_options.use_random_order = false;
  default_options.maximum_number_to_process = 0;
  default_options.number_computable_to_skip = 0;
  default_options.should_prefilter_files = false;
  % default_options.always_run_finished_compute_function = false;
  % default_options.should_filter_skipped = true;
  % default_options.use_all_file_ids = false;
  % default_options.run_all_file_ids = false;
  
  % Functions determining upon what we compute what:
  % This determines the files to load given a directory:
  % default_options.file_listing_function = @(given_load_location)dirr(given_load_location, '*.*');
  % default_options.file_listing_function = @(given_load_location)dirr([given_load_location, filesep, '*.stk']);
  function [listing] = default_file_listing_function(given_load_location)
    % [~, ~, listing] = dirr([given_load_location, filesep, '*.*'], 'name');
    [~, ~, listing] = dirr(given_load_location, 'name', '.+\..+');
  end
  default_options.file_listing_function = @(given_load_location)default_file_listing_function(given_load_location);
  % This interprets the files' paths (which might otherwise be done multiple times by should_load_function, should_compute_function, compute_function, and finished_compute_function):
  % default_options.file_id_function = @(given_file_listing){given_file_listing.name}.';
  default_options.file_id_function = @(given_file_listing)given_file_listing;
  % This determines whether this file id is of interest:
  default_options.should_load_function = @(given_file_id)true;
  % This determines how the file is converted to a Matlab variable:
  default_options.load_function = @(given_filename, given_file_id)load(given_filename);
  % This determines whether this file's data is of interest:
  default_options.should_compute_function = @(given_compute_function_arguments)true;
  % This computes data to be saved:
  default_options.compute_function = @(given_compute_function_arguments)[];
  % This postproceses the saved data, e.g., for visualization:
  default_options.finished_compute_function = @(given_compute_function_arguments)[];
  % This determines the filenames with which data are saved:
  default_options.save_filename_function = @(given_filenames, given_load_location, given_save_location)strrep(given_filenames, given_load_location, given_save_location);

  if ~exist('options', 'var')
    options = default_options;
  else
    [options] = process_options_structure(default_options, options);
  end

  
  filenames = options.file_listing_function(options.load_location)';
  filenames = regexprep(filenames, '/+', '/');
  file_ids = options.file_id_function(filenames);
  should_load_files = arrayfun(options.should_load_function, file_ids);
  D = false(size(should_load_files));

  features = [];
  % should_compute_files = false(length(filenames), 1);
  should_compute_files = true(length(filenames), 1);
  for file_index = 1:length(filenames)
    file_id = file_ids(file_index); 
    if should_load_files(file_index)
      continue
    end
    % should_compute_files(file_index) = true; 
    % This is often slow, but allow it with an option:
    if options.should_prefilter_files
      % current_load_location = [options.load_location, filenames{file_index}];
      current_load_location = filenames{file_index};
      % loaded_structure = load([options.load_location, filenames{file_index}]);
      loaded_structure = options.load_function(current_load_location, file_id);
      compute_function_arguments = struct(...
        'input_structure', loaded_structure ...
        , 'file_id', file_id ...
        , 'filename', current_load_location ...
        );
      if ~options.should_compute_function(loaded_structure)
        should_compute_files(file_index) = false; 
      end
    end
  end
  
  % if options.should_filter_skipped
    % filenames = filenames(~skipped_file_indices);
    % file_ids = file_ids(~skipped_file_indices); 
  % end
  % % filenames
  % skipped_file_indices = false(length(filenames), 1); 
  % % skipped_file_indices = true(length(filenames), 1); 
    
  
  % save_filenames = strrep(filenames, options.load_location, options.save_location);
  save_filenames = options.save_filename_function(filenames, options.load_location, options.save_location);
  
  % keyboard
  
  result = cell(length(filenames), 1);
  finished_values = cell(length(filenames), 1);
  
  file_indices = 1:length(filenames);
  if options.use_random_order
    file_indices = file_indices(randperm(length(file_indices)));
  end
  
  number_processed = 0;
  number_computable_skipped = 0;

  for file_index = file_indices
    file_id = file_ids(file_index); 
    % save_filename = filenames{file_index};
    % save_filename = save_filename(1:end - 4);
    save_filename = save_filenames{file_index};
    [base, name, extension] = fileparts(save_filename);
    % save_filename = [base, name];
    save_filename = name;
    mkdir(base)
    current_image_file = [options.load_location, filesep, save_filename];
    % current_save_location = options.save_location;
    current_save_location = base;
    % mkdir(current_save_location)
    
    % current_load_location = [options.load_location, filenames{file_index}];
    current_load_location = filenames{file_index};
    % fprintf(['current_load_location = ', current_load_location, '\n'])
    original_exists = exist(current_load_location, 'file');
    % if ~original_exists && ~options.run_all_file_ids
    if ~original_exists
      % skipped_file_indices(file_index)=true;
      continue
    end
    
    
    if options.dry_run
      can_start = true;
      final_exists = false;
    else
      % [can_start, final_name, final_exists] = chunk_start(current_save_location, save_filename, extension);
      [can_start, final_name, final_exists] = chunk_start(current_save_location, save_filename);
    end
    if ~can_start
      if ~final_exists
        if options.verbose, fprintf('%%%%%%%% Skipping ''%s''\n', save_filename), end
        continue
      else
        if options.verbose, fprintf('%%%%%%%% Loading result for ''%s''\n', save_filename), end
        current_result = load(final_name);
        % current_result = options.load_function(final_name, file_id);

      end
    end
    if options.verbose, fprintf('%%%%%%%% Can start ''%s''\n', save_filename), end
    processed_iteration = true;
    if options.number_computable_to_skip > number_computable_skipped
      number_computable_skipped = number_computable_skipped + 1;
      processed_iteration = false;
    else
      % skipped_file_indices(file_index)=false; 
      loaded_structure = options.load_function(current_load_location, file_id);

      compute_function_arguments = struct(...
        'input_structure', loaded_structure ...
        , 'file_id', file_id ...
        , 'filename', current_load_location ...
        );
      if ~options.should_compute_function(compute_function_arguments)
        should_compute_files(file_index) = false; 
      end
      if ~final_exists && ~options.just_run_finished_compute_function && ~isempty(options.compute_function) && should_compute_files(file_index)
        if options.verbose, fprintf('%%%%%%%% Computing ''%s''\n', save_filename), end
        % file_features = compute_function(loaded_segmentation_structure, file_id);
        current_result = options.compute_function(compute_function_arguments);
        % current_result
        if ~options.dry_run && ~options.load_only
          % save([current_save_location, filesep, save_filename], 'current_result')
          % save(final_name, 'current_result')
          % save(final_name, '-struct', 'current_result')
          save(final_name, 'current_result')
        end
        final_exists = true;
      end
    end
    
  
    if ~options.dry_run
      chunk_finish(current_save_location, save_filename);
    end

    % if final_exists || options.always_run_finished_compute_function
    % if final_exists
    if final_exists && processed_iteration
      if options.verbose, fprintf('%%%%%%%%%% Calling options.finished_compute_function!\n'), end
      compute_function_arguments.result = current_result;
      % options.finished_compute_function(compute_function_arguments);
      current_finished_result = options.finished_compute_function(compute_function_arguments);
      % current_finished_result
      finished_values{file_index} = current_finished_result;
    end
    
    % if ~options.just_save_results
    % if ~options.just_save_results && ~options.just_run_finished_compute_function
    % if ~options.just_save_results && final_exists
    if ~options.just_save_results && final_exists && processed_iteration
      result{file_index} = current_result;
    end
    
    % number_processed = number_processed + 1;
    % if ~skipped_file_indices(file_index)
    if processed_iteration
      number_processed = number_processed + 1;
    end
    % processed_iteration
    % number_processed
    if options.maximum_number_to_process > 0 && number_processed >= options.maximum_number_to_process
      break
    end
  end
  %fprintf('done\n')
  %fprintf('\n')
  
  % if options.should_filter_skipped
    % filenames = filenames(~skipped_file_indices);
    % file_ids = file_ids(~skipped_file_indices); 
    % result = result(~skipped_file_indices, :); 
    % finished_values = finished_values(~skipped_file_indices, :); 
  % end

  values = result;
  whos values file_ids finished_values
  
  
end
