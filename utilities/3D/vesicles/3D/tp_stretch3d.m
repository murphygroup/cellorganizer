function im2 = tp_stretch3d(img, zsize, method)
% TP_STRETCH3D resize a 3D stack along the z-direction using trilinear
% interpolation

% Author: Tao Peng
% Edited: Ivan E. Cao-Berg 
% 6/11/13 D. Sullivan added support for multiple interpolation methods
%
% Copyright (C) 2011-2012 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
% Carnegie Mellon University
%
% July 26, 2012 Devin S. Changed method to perform bilinear interpolation at resize
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

%D. Sullivan 6/11/13 added method choice
if nargin<3
    method = 'bilinear';
elseif isempty(method)
    method = 'bilinear';
end

if size(img,3) == zsize
    im2 = img;
    return;
end

im2 = zeros(size(img,1),size(img,2),zsize);
for i = 1:size(img,1)
    I = squeeze(img(i,:,:));
    %devins 26/7/2012
    %I = imresize(I,[size(I,1) zsize]);
    %D. Sullivan 6/11/13 added method choice
%     I = imresize(I,[size(I,1) zsize],'bilinear');
    I = imresize(I,[size(I,1) zsize],method);
    im2(i,:,:) = I;
end
