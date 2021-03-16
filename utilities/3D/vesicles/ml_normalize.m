function NormalImage = ml_normalize(img);

% ML_NORMALIZE normalizes an image by adjusting its major angle and center
%   NORMALIMAGE=ML_NORMALIZE(IMAGE) return an image with major angle
%   along image rows and weight center at the center of the image.
%   The size of image will also be changed to the next power of 2.
%   For example, an image with size 500x200 will return a 512x512 image.

% Copyright (C) 2006  Murphy Lab
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

% Written by Adrienne Wells
% Jun. 29, 2005 Modified by T. Zhao

Theta=ml_majorangle(img)*180/pi;

img=ml_rotate(img,-Theta);

NormalImage = ml_translate(img);
