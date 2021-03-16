function img2 = ml_imcropbg(img,threshold)
%ML_IMCROPBG Crop part of background out
%   IMG2 = ML_IMCROPBG(IMG) returns an image by cropping background of IMG.
%   Here background means pixels with value no greater than 0.
%   
%   IMG2 = ML_IMCROPBG(IMG,THRESHOLD) defines background by a specified
%   value THRESHOLD.

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

%   20-Oct-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('1 or 2 arguments are required')
end

if nargin < 2
    threshold = 0;
end

img2=img;
if length(size(img))>2
    img = sum(img,3);
end
img(img>threshold) = 1;
img(img<=threshold) = 0;

colidx = find(sum(img,1)>0);
rowidx = find(sum(img,2)>0);

if ~isempty(colidx)
    if length(colidx)==1
        img2 = img2(:,colidx,:);
    else
        img2(:,[1:colidx(1)-1 colidx(end)+1:end],:) = [];
    end
end

if ~isempty(rowidx)
    if length(rowidx)==1
        img2 = img2(rowidx,:,:);
    else
        img2([1:rowidx(1)-1 rowidx(end)+1:end],:,:) = [];
    end
end

% img2(:,sum(img,1)==0)=[];
% img2(sum(img,2)==0,:)=[];
