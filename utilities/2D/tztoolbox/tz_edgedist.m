function img2=tz_edgedist(img)
%TZ_EDGEDIST Build a distance map for the edge of an image.
%   IMG2 = TZ_EDGEDIST(IMG) returns a distance map of the edge in IMG, of
%   which pixels with intensity 0 are taken as background.
%   
%   See also

%   ??-???-???? Initial write TINGZ
%   31-OCT-2004 Modified TINGZ
%       - add comments
%   17-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

% Copyright (C) 2007  Murphy Lab
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

img2=img>0;
img2=bwperim(img2);
img2=bwdist(img2);
