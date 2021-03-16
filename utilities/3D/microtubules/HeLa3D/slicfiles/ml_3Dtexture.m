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

function textf = xc_3Dtexture(img)
% FUNCTION TEXTF = XC_3DTEXTURE(IMG)
% Calculate 3D version texture features.  The major difference is that 
% gray-level cooccurence matrices are build on 13 directions (instead of 4 
% in 2D images).
%
% img: input 3D image.  Background subtraction recommended.
% textf: a 14 * 15 matrix texture features; Rows 1 - 14 are 14 statistics 
% defined by Haralick and columns 1 - 13 are 13 different direction.  The 
% 14th column is the average of the 13 directions. Column 15 is the range
% across the 13 directions.
