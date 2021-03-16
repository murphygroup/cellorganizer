function branch_points = ml_find_branch_points( img)
% function branch_points = ml_find_branch_points( img)
% img is the image with the skeleton of one object
% bran_points is the image with the branch points of the skeleton of that object

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

% Written by Meel Velliste

kernel = [0 1 0; 1 0 1; 0 1 0];
branch_points = conv2(img,kernel,'same').*img;
branch_points(find(branch_points<3))=0;
branch_points(find(branch_points))=1;
