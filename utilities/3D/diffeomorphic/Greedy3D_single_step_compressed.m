function [d_source, d_target, d_source_distance, d_target_distance] = ...
  Greedy3D_single_step_compressed(compressed_source, compressed_target, options)
  % Internal function for Greedy3D_lambda_pre_compressed.
  % 
  % Computes the derivatives of the deformation fields and distances for the (currently deformed) source and target images.
  %
  % 2013-02-17 tebuck: Changing the default options to those that are currently most useful or are being used in CellOrganizer.
  
  % Default option values and options processing (these are now set in Greedy3D_compressed_default_options.m):

  % default_options.single_sided = false; 
  % default_options.periodic_space = false; 
  % default_options.filter_radius = 32; 
  % default_options.verbose = 1; 
  default_options = Greedy3D_compressed_default_options();
  if ~exist('options', 'var')
    options = default_options;
  else
    options = process_options_structure(default_options, options);
  end

  % d_source = repmat({cell(size(compressed_source))}, [3, 1]); 
  d_source = repmat({compressed_source}, [3, 1]); 
  d_target = []; 
  if (~options.single_sided)
    d_target = d_source; 
  end
  
  window_mode = 'circular'; 
  if (~options.periodic_space)
    window_mode = 'replicate'; 
  end
  
  
  d_source_distance = 0; 
  d_target_distance = 0; 
  % For each window, pad it with paddarray in replicate mode with
  % filter_radius + 1 elements in the first two dimensions and
  % replace the padding with portions of the adjacent windows, if
  % possible. Run single_step_window on the padded window, store
  % its returned derivatives in compressed form, and compute the
  % distance moved by the deformation. Note that this assumes that
  % windows are larger than filter_radius + 1.
  [number_rows, number_columns] = size(compressed_source.image); 
  perform_offset = true; 
  z_size = []; 
  offsets = []; 
  for row = 1:number_rows
    for column = 1:number_columns
      if options.verbose > 1
        fprintf([mfilename ': r%05d,c%05d\n'], row, column);
      end
      % Pad with an extra pixel because of del2:
      %whos
      padded_window_source = compressed_source.get_padded_window(...
        row, column, ceil(options.filter_radius + 1), window_mode);
      padded_window_target = compressed_target.get_padded_window(...
        row, column, ceil(options.filter_radius + 1), window_mode);
      z_size = size(padded_window_source, 3); 
      % del2 consumes an extra pixel around the edges, so we pass
      % the actual radius of the filter then are returned only the
      % valid area (filter_radius + 1) smaller on all sides in XY:
      [ds, dt, dsd, dtd, offsets] = ...
        Greedy3D_single_step_window(padded_window_source, padded_window_target, offsets, options);
      for ind = 1:3
        d_source{ind} = d_source{ind}.set_window(row, column, ds{ind});
        if (~options.single_sided)
          d_target{ind} = d_target{ind}.set_window(row, column, dt{ind});
        end
      end
    end
  end
  % Keep velocity zero at the edges:
  for dim_ind = 1:3
    d_source{dim_ind} = ...
      d_source{dim_ind}.set_borders(0, z_size == 1);
    if (~options.single_sided)
      d_target{dim_ind} = ...
        d_target{dim_ind}.set_borders(0, z_size == 1);
    end
  end
  d_source_distance = WindowedImage.LDDMM_distance(d_source);
  if (~options.single_sided)
    d_target_distance = WindowedImage.LDDMM_distance(d_target);
  end

 