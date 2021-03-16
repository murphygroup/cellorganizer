function [frames] = generate_frame_from_shape_space_model(model, options)
  % Generate a frame from a trained shape space model. Return value is a structure with fields interpolated_image, total_wall_time, and total_cpu_time.
  % 
  % % Dependencies from File Exchange:
  % % inhull
  % % sort_nat
  % % export_fig
  % % CompressionLib
  % % convnfft
  % 
  % Dependencies from File Exchange:
  % DataHash
  % Dependencies from the Murphy and Rohde labs:
  % peng_interpolation_code
  % 
  % 2012-09-12 tebuck: Copied from test_walk_ternary023c.m.
  % 2013-07-02 tebuck: Modifying to take functions instead of assuming filenames.
  % 2013-07-26 gj:     Modified to work with cell organizer
  
  shape_space = model.nuclearShapeModel;

  
%   options = process_options_structure(shape_space.shape_space_options, options); 
  positions = options.positions;
  options = model.nuclearShapeModel.shape_space_options;
  options.positions = positions;
  options.image_function = model.nuclearShapeModel.imfunc;
  options.use_known_distances = false;
  
  % options, error
  
  registration_options = options;
  
  
  shape_space2 = {shape_space.positions, shape_space.convex_hull, shape_space.tessellation};
  if options.use_known_distances
    shape_space2{end + 1} = shape_space.distances;
  end
  frames = cell(size(options.positions, 1), 1);
  for position_index = 1:length(frames)
    frames{position_index} = render_point_3d_windowed(options.positions(position_index, :), shape_space2, options.image_function, 0, registration_options);
  end


end

