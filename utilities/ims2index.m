function indexed = ims2index( img, mapping )
%IMG2INDEX Converts an intensity image into an indexed image.
%
%Input arguments    Description
%img                A cell array of multidimensional images of the same size
%mapping        	Indices to assign to each image in img (to control 
%                   precedence of one pattern over another)

% Author: Michelle Mackie (mmackie@andrew.cmu.edu)
% June 6, 2012 M. Mackie Removed parameter 'indexOption', only uses
%                        'priority'
% June 19, 2012 I. Cao-Berg Removed top header and added a check of the
%                           input arguments
% July 26, 2012 R.F. Murphy Fix output for last image; add mapping to
%                           control priority of indices
%
% Copyright (C) 2012 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
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

%icaoberg june 19, 2012
indexed = [];
if isempty( img )
    return
end

if nargin < 2
    mapping = eye(length(img),1)
end

for i=1:length(img)
    image = double(img{i});
    image(image>0) = mapping(i);
    %check to see if sizes are equal
    if ~isequal(size(image),size(img{end}))
        warning('image sizes are not equal');
        return;
    end
    img{i} = image;
end

indexed = zeros(size(img{i}));


for i=1:length(img)
    %check to see if images are empty
    if isempty(img{i})
        warning('image is empty');
        return;
    end
    indexed = indexed+img{i};
    indexed(indexed>i)=i;
end
% end
