function img = im2reshape( img )
% IM2RESHAPE Reshapes a 3D images into a 2D representation
%
% List Of Input Arguments     Descriptions
% -----------------------     ------------
% img                         a 3D binary or realvalued image
%
% List Of Outputs     Descriptions
% ---------------     ------------
% img                 a 2D representation of the image

% Author: Ivan Cao-Berg
%
% Copyright (C) 2013-2016 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
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

if size( img, 3 ) > 1
    img = reshape( img, size( img, 1 ), [] );
else
    warning( 'Image array must 3D' );
    img = [];
end