function [t_cell_info] = tcell_get_template_info(t_cell_info, options)
  % 2016-03-06 xruan: Copied from master_script_tcell_get_template_info.m.
  % 
  % 
  % 
  % 
  % 
  % % Dependencies:
  % % From the File Exchange: export_fig, pmkmp, toolbox_graph, toolbox_fast_marching, affine, Snake3D, Mesh_voxelisation, dirr, inhull
  % % From tebuck: modified version of Jieyue Li's HPA_lib, rasterize_mesh

  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Figure out how large the window surrounding a synapse point should be and create a template cell shape of similar size:

  
  % Variables referenced more than once can be copied to local variables:
  cell_diameter = t_cell_info.options.cell_diameter;

  
  room = round(cell_diameter * t_cell_info.options.room_scale);
  % Z size will be set to the Z size of each original image:
  window_size_2d = (room * 2 + 1) .* [1, 1];
  window_center_2d = (room + 1) .* [1, 1];
  % maximum_number_slices = 0;
  % template_number_slices = ceil(room * 1.25);
  template_number_slices = room;
  
  default_template_options = get_default_template_options();
  template_options = default_template_options;
  template_options.imx = window_size_2d(2);
  % template_options.imx = room;
  template_options.imy = template_options.imx;
  template_options.imz = template_number_slices;
  template_options.xc = template_options.imx / 2 + .5;
  template_options.yc = template_options.imy / 2 + .5;
  template_options.zc = template_options.imz / 2 + .5;
  template_options.xr = cell_diameter;
  template_options.yr = cell_diameter / 2;
  template_options.zr = cell_diameter / 2;
  % template_options
  
  template_image = get_template(template_options);
  template_cell_radius = (template_options.yr + template_options.zr) / 2;
  template_synapse = [template_options.yc, template_options.xc, template_options.zc];
  % template_centroid = template_synapse + [0, 3 / 16 * template_options.xr, 0]
  empirical_template_centroid = image_centroid(template_image);
  template_centroid = empirical_template_centroid;
  % XYZ order, not RCZ, so _dimensions instead of _size:
  template_dimensions = [template_options.imy, template_options.imx, template_options.imz];
  template_size = template_dimensions([2, 1, 3]);
  
  % if t_cell_info.options.use_two_point_synapses
    % Produce two-point landmarks:
    template_synapse_segment = repmat([template_options.yc, template_options.xc, template_options.zc], [2, 1]);
    template_synapse_segment(:, 1) = template_synapse_segment(:, 1) + [-1; 1] .* template_options.yr;
  % end

  % This "cropping" might just be a variable naming relic:
  cropped_size = template_size;
  
  % This "cropping" is only used after the models have been built:
  template_crop_x = sum(sum(template_image, 1), 3) > 0;
  template_crop_y = sum(sum(template_image, 2), 3) > 0;
  template_crop_z = sum(sum(template_image, 1), 2) > 0;
  template_crop_function = @(given_image)given_image(template_crop_y, template_crop_x, template_crop_z);
  
  
  function [result_background_level] = background_level_function(given_image)
    result_background_level = mode(given_image(given_image < mean(given_image(:))));
  end
  
  
  function [result_image] = background_subtraction_function(given_image)
    background_level = background_level_function(given_image);
    result_image = given_image - background_level;
    result_image(result_image < 0) = 0;
  end
  
  
  % Is this BK's X, where X is row?
  template_occupied_x_slices = sum(sum(template_image, 2), 3) > 0;
  template_occupied_x_slice_pad_function = @(x)padarray(padarray(reshape(x, [], 1), [find(template_occupied_x_slices, 1, 'first') - 1, 0, 0], 'pre'), [length(template_occupied_x_slices) - find(template_occupied_x_slices, 1, 'last'), 0, 0], 'post');
  
  template_occupied_y_slices = sum(sum(template_image, 1), 3) > 0;
  template_occupied_y_slice_pad_function = @(x)padarray(padarray(reshape(x, [], 1), [0, find(template_occupied_y_slices, 1, 'first') - 1, 0], 'pre'), [0, length(template_occupied_y_slices) - find(template_occupied_y_slices, 1, 'last'), 0], 'post');
  
  approximate_synapse_slice = round(prctile(find(template_occupied_x_slices), 10));
  approximate_synapse_slice_template_image = false(size(template_image));
  approximate_synapse_slice_template_image(approximate_synapse_slice, :, :) = true;
  approximate_synapse_slice_template_image = approximate_synapse_slice_template_image & template_image;
  
  % 2013-10-03: For current template parameters, Full Stimulus_ARP3 shows about 6/23 slices, Full Stimulus_Actin 5-6, Full Stimulus_CPalpha1 5-7. Let's say up to 25% to be conservative and mostly exclude the nucleus:
  approximate_cytoplasm_slices = round(prctile(find(template_occupied_x_slices), [0, 25]));
  approximate_cytoplasm_slices = approximate_cytoplasm_slices(1):approximate_cytoplasm_slices(2);
  approximate_cytoplasm_slice_template_image = false(size(template_image));
  approximate_cytoplasm_slice_template_image(approximate_cytoplasm_slices, :, :) = true;
  approximate_cytoplasm_slice_template_image = approximate_cytoplasm_slice_template_image & template_image;
  
  % 2014-06-06: Get coordinates for positions in the cell in BK's order:
  [template_x, template_y, template_z] = ndgrid(1:size(template_image, 1), 1:size(template_image, 2), 1:size(template_image, 3));
  % Normalized by template radius:
  template_cylindrical_r = sqrt(((template_y - template_options.yc) ./ template_options.yr).^2 + ((template_z - template_options.zc) ./ template_options.zr).^2);
  template_cylindrical_theta = atan2(template_z - template_options.zc, template_y - template_options.yc);
  
  
  t_cell_info.preprocessing_info = struct();
  t_cell_info.preprocessing_info.room = room;
  t_cell_info.preprocessing_info.window_size_2d = window_size_2d;
  t_cell_info.preprocessing_info.window_center_2d = window_center_2d;
  t_cell_info.preprocessing_info.cropped_size = cropped_size;
  t_cell_info.preprocessing_info.background_level_function = @background_level_function;
  t_cell_info.preprocessing_info.background_subtraction_function = @background_subtraction_function;
  
  t_cell_info.template_info = struct();
  t_cell_info.template_info.template_options = template_options;
  t_cell_info.template_info.template_number_slices = template_number_slices;
  t_cell_info.template_info.template_image = template_image;
  t_cell_info.template_info.template_cell_radius = template_cell_radius;
  t_cell_info.template_info.template_synapse = template_synapse;
  t_cell_info.template_info.empirical_template_centroid = empirical_template_centroid;
  t_cell_info.template_info.template_centroid = template_centroid;
  t_cell_info.template_info.template_dimensions = template_dimensions;
  t_cell_info.template_info.template_size = template_size;
  t_cell_info.template_info.template_synapse_segment = template_synapse_segment;
  t_cell_info.template_info.template_crop_x = template_crop_x;
  t_cell_info.template_info.template_crop_y = template_crop_y;
  t_cell_info.template_info.template_crop_z = template_crop_z;
  t_cell_info.template_info.template_crop_function = template_crop_function;
  
  t_cell_info.template_info.template_occupied_x_slices = template_occupied_x_slices;
  t_cell_info.template_info.template_occupied_x_slice_pad_function = template_occupied_x_slice_pad_function;
  t_cell_info.template_info.template_occupied_y_slices = template_occupied_y_slices;
  t_cell_info.template_info.template_occupied_y_slice_pad_function = template_occupied_y_slice_pad_function;
  t_cell_info.template_info.approximate_synapse_slice = approximate_synapse_slice;
  t_cell_info.template_info.approximate_synapse_slice_template_image = approximate_synapse_slice_template_image;
  t_cell_info.template_info.approximate_cytoplasm_slices = approximate_cytoplasm_slices;
  t_cell_info.template_info.approximate_cytoplasm_slice_template_image = approximate_cytoplasm_slice_template_image;
  t_cell_info.template_info.template_x = template_x;
  t_cell_info.template_info.template_y = template_y;
  t_cell_info.template_info.template_z = template_z;
  t_cell_info.template_info.template_cylindrical_r = template_cylindrical_r;
  t_cell_info.template_info.template_cylindrical_theta = template_cylindrical_theta;
  
  
end  

