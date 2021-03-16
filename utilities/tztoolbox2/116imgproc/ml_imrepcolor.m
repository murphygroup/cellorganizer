function img2 = ml_imrepcolor(img,color1,color2)
%ML_IMREPCOLOR Replace color in a rgb image.
%   IMG2 = ML_IMREPCOLOR(IMG) replace the color in the rgb image IMG from
%   COLOR1 to COLOR2. 
%   
%   See also

%   07-Nov-2006 Initial write T. Zhao
%   Copyright (c) 2006 Murphy Lab
%   Carnegie Mellon University
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation; either version 2 of the License,
%   or (at your option) any later version.
%   
%   This program is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   General Public License for more details.
%   
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
%   02110-1301, USA.
%   
%   For additional information visit http://murphylab.web.cmu.edu or
%   send email to murphy@cmu.edu


if nargin < 1
    error('Exactly 1 argument is required');
end

idx = find(img(:,:,1)==color1(1) & img(:,:,2)==color1(2) & ...
           img(:,:,3)==color1(3));

for i=1:3
    img2 = img(:,:,i);
    img2(idx) = color2(i);
    img(:,:,i) = img2;
end

img2 = img;
