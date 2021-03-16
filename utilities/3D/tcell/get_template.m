function [templ, cross_section_radii] = get_template(options)
  % Generate a binary elipsoid as a mask for morphing.
  % 
  % Created 2011-11-29 tebuck.

  % Process options:
  default_options = get_default_template_options(); 
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
  %options

  imx = options.imx; 
  imy = options.imy; 
  imz = options.imz; 
  xc = options.xc; 
  yc = options.yc; 
  zc = options.zc; 
  xr = options.xr; 
  yr = options.yr; 
  zr = options.zr; 
  
  % From Baek Hwan Cho 2011-11-29:
  %templ = zeros(imx,imy, imz);
  templ = false(imx,imy, imz);
  for x=1:imx
     for y=1:imy
         for z=1:imz
             if (x-xc)*(x-xc)/(xr*xr)+(y-yc)*(y-yc)/(yr*yr)+(z-zc)*(z-zc)/(zr*zr)<=1 && x>=xc
                 %templ(x,y,z)=1;
                 templ(x,y,z)=true;
             end
         end
     end
  end

  if nargout > 1
    cross_section_radii = nan(imx, 2); 
    minimum_x = ceil(xc); 
    maximum_x = floor(xc + xr); 
    if maximum_x - xc == xr
      maximum_x = maximum_x - 1;
    end
    x_subtractions = ((minimum_x:maximum_x)-xc).*((minimum_x:maximum_x)-xc)./(xr*xr); 
    % Ask what maximum y - yc gives <= 1 in template generation
    % loop at this x and z = zc (x is row, y is column here):
    cross_section_radii(minimum_x:maximum_x, 1) = yr .* sqrt(1 - x_subtractions); 
    % Then for z:
    cross_section_radii(minimum_x:maximum_x, 2) = zr .* sqrt(1 - x_subtractions); 
    %templ
    %cross_section_radii
  end
  %plot(cross_section_radii), pause
