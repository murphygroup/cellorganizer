function res = OME_getResolution(img_path)
%OME_getResolution returns the resolution found in img_path's metadata
%
%Inputs:
% img_path = string
%
%Outputs:
% res = 1xN vector of [x y] or [x y z]

%Author: Tim Majarian 7/28/16
% Copyright (C) 2016  Murphy Lab
% Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.
%
% For additional information visit http://murphylab.web.cmu.edu or
% send email to murphy@cmu.edu

info = bfopen(img_path);
info = info{1,4};
xres = double(info.getPixelsPhysicalSizeX(0).value());
yres = double(info.getPixelsPhysicalSizeY(0).value());
zres = double(info.getPixelsPhysicalSizeZ(0).value());
if exist( 'zres', 'var' )
    res = [xres yres zres];
else
    res = [xres yres];
end
end