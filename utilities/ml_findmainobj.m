function [ imobj ] = ml_findmainobj(img, conn )
%finds the largest non-zero pixel region in an image

% grj 3/29/13
% Copyright (C) 2007-2013  Murphy Lab
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

if ~exist('conn', 'var')
    if ndims(img) == 2
        conn = 8;
    else
        conn = 26;
    end
end

disp('Finding objects')
tic
objs = ml_findobjs(img, conn);
toc

[~, ind] = max(cellfun(@(x) size(x, 1), objs));

if ~isempty( ind )
    vox = objs{ind};
else
    warning([num2str(length(unique(img))) ' unique values found in image.' ])
    warning( 'No objects found in image. Returning input image.' );
    imobj = img;
    return
end

pix = vox(:,4);
vox = double(vox(:,1:3));

inds = sub2ind(size(img), vox(:,1), vox(:,2), vox(:,3));

imobj = zeros(size(img));
imobj(inds) = pix;

end

