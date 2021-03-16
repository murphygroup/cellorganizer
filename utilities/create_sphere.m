function s = create_sphere(radius,x,y,z)
%MYSPHERE Helper method that creates a 3d sphere image
%
%Input: size of the image (x,y,z), and radius of the sphere(radius > 10)
%Output: 3d image

% Author: Yue Yu(yuey1@andrew.cmu.edu)
%
% Copyright (C) 2012 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
% Carnegie Mellon University
%
% October 12, 2012 I. Cao-Berg Documented method, added check of input
% arguments, encapsulated method in try/catch method so that method returns
% empty array if it fails
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
[X,Y,Z] = ndgrid(-x/2:x/2,-y/2:y/2,-z/2:z/2);
 sphere = (X.^2+Y.^2+Z.^2 <= radius^2);
 sp = double(sphere);
 s = double(bwperim(sp));
 %trim the top and bottom of the sphere
 s(:,:,1:z/2-radius+5) = 0;
 s(:,:,z/2+radius-5:end) = 0;
 
                
                
