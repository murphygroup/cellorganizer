function img2 = ml_imcropbg3d(img,threshold)
%ML_IMCROPBG3D Crop background of a 3D image.
%   IMG2 = ML_IMCROPBG3D(IMG) returns a smaller 3D image of the input 3D
%   image IMG. In the returned image the background is croped. Here
%   background means pixels with value no greater than 0.
%   
%   IMG2 = ML_IMCROPBG3D(IMG,THRESHOLD) defines background by a specified
%   value THRESHOLD.
%
%   See also

%   11-Oct-2006 Initial write T. Zhao
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

if nargin < 2
    threshold = 0;
end

img2 = img;
img(img>threshold) = 1;
img(img<=threshold) = 0;

img2(sum(sum(img,3),2)==0,:,:)=[];
img2(:,sum(sum(img,3),1)==0,:)=[];
img2(:,:,sum(sum(img,2),1)==0)=[];
