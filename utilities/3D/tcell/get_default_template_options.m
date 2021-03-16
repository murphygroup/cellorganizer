function [default_options] = get_default_template_options()
  % Baek-Hwan Cho's 
  % 
  % Created 2011-12-19 tebuck.

  % Process options:
  default_options = struct(); 
  default_options.imx = 50;    % size of template image
  default_options.imy = 50;
  default_options.imz = 46;
  default_options.xc = 10;      % center coordinates and size of half elipsoid(template)
  default_options.yc = 25.5;
  default_options.zc = 23.5;
  default_options.xr = 22;
  default_options.yr = 11;
  default_options.zr = 11;
