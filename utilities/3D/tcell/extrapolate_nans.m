function [result_image] = extrapolate_nans(given_image, options)
  % Do a distance transform for all NaN values and fill them in with their nearest neighbors.
  %
  % Dependencies:
  % File Exchange: http://www.mathworks.com/matlabcentral/fileexchange/4551-inpaintnans
  %
  % Tests:
  % a = [1, 2, 3, nan, nan, 6, nan, 8, nan], extrapolate_nans(a, struct('method', 'diffusion', 'diffusion_method', 0))
  % a = [1, 2, 3, nan(1, 10), 6, 7, nan, 9, 10, nan], ae = extrapolate_nans(a, struct('method', 'diffusion', 'diffusion_method', 1)), plot(a, 'r-'), hold on, plot(ae, 'b:'), hold off
  % a = [1, 2, 3, nan(1, 10), 6, 7, nan, 9, 10, nan], ae = extrapolate_nans(a, struct('method', 'diffusion', 'diffusion_method', 4)), plot(a, 'b-'), hold on, plot(ae, 'r:'), hold off
  %
  % 2013-05-01 tebuck: Created.
  % 2013-07-09 tebuck: Added options.method = 'diffusion', uses FEX 4551. Should eventually use FEX 34356 after making parameterized functions for it.
  % 2013-09-28 tebuck: Added options.diffusion_method to set the method used by FEX 4551.
  
  default_options = struct();
  default_options.method = 'nearest';
  default_options.diffusion_method = 2;

  if ~exist('options', 'var')
    options = default_options; 
  else
    options = process_options_structure(default_options, options);
  end
  % options

  switch options.method
  
    case 'diffusion'
      result_image = inpaint_nans(given_image, options.diffusion_method);
    
    case 'nearest'
      result_image = given_image;
      bad_pixels = isnan(given_image);
      if all(bad_pixels(:))
        return
      end
      while any(bad_pixels(:))
        [bad_pixel_distances, neighbors] = bwdist(~bad_pixels);
        % adjacent_bad_pixels = bad_pixel_distances == 1;
        adjacent_bad_pixels = bad_pixels;
        result_image(adjacent_bad_pixels) = result_image(neighbors(adjacent_bad_pixels));
        bad_pixels(adjacent_bad_pixels) = false;
      end
      
  end
end