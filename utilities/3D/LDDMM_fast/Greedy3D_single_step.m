function [d_source, d_target, d_source_distance, d_target_distance] = ...
  Greedy3D_single_step(source, target, options)
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
    options = process_options_structure_fast(default_options, options);
  end

  % d_source = repmat({cell(size(compressed_source))}, [3, 1]); 
  d_source = repmat({source}, [3, 1]); 
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
  [number_rows, number_columns] = size(1); 
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
      padded_window_source = get_padded_window(...
        source, row, column, ceil(options.filter_radius + 1), window_mode);
      padded_window_target = get_padded_window(...
        target, row, column, ceil(options.filter_radius + 1), window_mode);
      z_size = size(padded_window_source, 3); 
      % del2 consumes an extra pixel around the edges, so we pass
      % the actual radius of the filter then are returned only the
      % valid area (filter_radius + 1) smaller on all sides in XY:
      [ds, dt, dsd, dtd, offsets] = ...
        Greedy3D_single_step_window_fast(padded_window_source, padded_window_target, offsets, options);
      for ind = 1:3
        d_source{ind} = ds{ind};
        if (~options.single_sided)
          d_target{ind} = dt{ind};
        end
      end
    end
  end
  % Keep velocity zero at the edges:
  for dim_ind = 1:3
    d_source{dim_ind} = ...
      set_borders(d_source{dim_ind}, 0, z_size == 1);
    if (~options.single_sided)
      d_target{dim_ind} = ...
        set_borders(d_target{dim_ind}, 0, z_size == 1);
    end
  end
  if z_size > 1
      d_source_distance = sqrt(mean(dsd .^ 2));
      % d_source_distance = LDDMM_distance(d_source)
      if (~options.single_sided)
        d_target_distance = sqrt(mean(dtd .^ 2));
      end
  else
      d_source_distance = dsd;
      d_target_distance = dtd;
  end
end
 