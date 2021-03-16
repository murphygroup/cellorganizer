function [slice_highest_orders] = get_template_highest_orders(highest_order, options)
  % Process options:
  default_options = get_default_template_options(); 
  default_options.highest_order_same_across_slices = false;
  if ~exist('options', 'var') || isempty(options)
    options = default_options;
  else
    option_names = fieldnames(default_options); 
    for index = 1:length(option_names)
      current_option = option_names{index}; 
      if ~isfield(options, current_option)
        options = setfield(options, current_option, getfield(default_options, current_option)); 
      end
    end
  end
  %options
  
  [templ, cross_section_radii] = get_template(options);
  % BHC's convention is row is x, column is y:
  occupied_slices = find(isfinite(cross_section_radii(:, 1)'));
  if options.highest_order_same_across_slices
    % slice_highest_orders = highest_order .* ones(length(cross_section_radii), 1);
    slice_highest_orders = highest_order .* ones(length(occupied_slices), 1);
  else
    % slice_highest_orders = ceil(highest_order .* min(cross_section_radii(occupied_slices, :) ./ repmat([options.yr, options.zr], length(occupied_slices), 1), [], 2));
    slice_highest_orders = round(highest_order .* min(cross_section_radii(occupied_slices, :) ./ repmat([options.yr, options.zr], length(occupied_slices), 1), [], 2));
    if false
      % Debug:
      slice_highest_orders
      pause
    end
  end
  