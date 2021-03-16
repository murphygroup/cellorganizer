function answer = has_rois( img_path )
% Input
% img (a valid OME.TIFF file)
% Output
% true if the number of ROIs > 0
%
% Xin Lu (xlu2@andrew.cmu.edu)
%
% Copyright (C) 2017-2018 Murphy Lab
% Computational Biology Department
% Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.o
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

number_of_rois = get_number_of_rois( img_path );

if number_of_rois <= 0
    answer=false;
else
    answer=true;
end
end