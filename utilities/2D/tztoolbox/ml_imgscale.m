
function scaledimage = ml_imgscale(image)
% ML_IMGSCALE scales the range of values in the input image to 0-1
%   ML_IMGSCALE(IMAGE) scales the values in IMAGE to be between 
%	0 and 1.
%   If you have the image processing toolbox, please use mat2gray 
%   instead of this function.

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

% 07 Aug 98
% Written by Michael Boland
% $Id: ml_imgscale.m,v 1.3 2006/06/27 13:33:47 tingz Exp $

if ~image
	error('Invalid input image') ;
end

%
% Subtract the minimum value and divide by the resulting maximum
%
scaledimage = image - min(image(:)) ;
scaledimage = scaledimage/max(scaledimage(:)) ;


