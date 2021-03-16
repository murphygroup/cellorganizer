function [result_image] = draw_text(lines, base_image, position, options)
% Requires bitmap_plot_version2 from the MathWorks File Exchange.
% 2012-03-16 tebuck: Created.
% 2012-11-10 tebuck: Added real-valued dilation option to make text more legible.

  % From bitmaptext.m:
  Font_sizes_x=[8 10 11];
  
  if ischar(lines)
    lines = {lines};
  end
  
  % Process options:
  default_options = struct(); 
  default_options.FontSize = 1;
  % Zero is horizontal, left to right, one is vertical bottom to
  % top, others unimplemented as yet:
  default_options.direction = 0;
  % default_options.centered = false;
  % -1 for left, 0 for centered, 1 for right:
  default_options.alignment = -1;
  default_options.dilation = 0;

  if ~exist('options', 'var')
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
  warning('off', 'register_images:unknownoption');
  % options

  if options.dilation > 0
    result_image = zeros(size(base_image)); 
  else
    result_image = base_image; 
  end

% $$$   if options.direction == 0
% $$$     position = fliplr(position); 
% $$$   end
  
  if options.direction == 1
    result_image = flipdim(permute(result_image, [2, 1, 3]), 2);
    position(2) = size(base_image, 1) - position(2); 
    position = fliplr(position); 
  end
  
  % Properly align text:
  % if options.centered
    % position(1) = position(1) - Font_sizes_x(options.FontSize) * .5 * max(cellfun(@(x)length(x), lines));
  % end
  position(1) = position(1) - Font_sizes_x(options.FontSize) * (options.alignment + 1) * .5 * max(cellfun(@(x)length(x), lines));
  
  position = round(position);
  
% $$$   if options.direction == 0
% $$$     position = fliplr(position); 
% $$$   end
  position = fliplr(position); 

  result_image = bitmaptext(lines, result_image, position, options);

  if options.dilation > 0
    % Approximately antialiased disc:
    ks = options.dilation + 1.; kw = ceil(ks + 1) * 2 + 1; [kx, ky] = meshgrid(1:kw, 1:kw); kc = ceil(kw / 2 + .5); k = min(max(ks - sqrt((kx - kc) .^ 2 + (ky - kc) .^ 2) ./ (1), 0), 1);
    result_image = imdilate(result_image, strel('arbitrary', k > 0, k)) - 1;
    % This should use options.Color eventually:
    result_image = base_image .* (1 - result_image) + result_image .* 1;
  end

  
  if options.direction == 1
    result_image = ipermute(flipdim(result_image, 2), [2, 1, 3]);
  end
  
  
