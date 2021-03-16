function [I] = convert_contour_to_image(cell_landmarks, nuc_landmarks, imageSize)
% convert contour of cell/nuc shape to image. 

% Author: Xiongtao Ruan (xruan@andrew.cmu.edu)
%
% Copyright (C) 2013-2017 Murphy Lab
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin < 3
    if ~isempty(cell_landmarks)
        imageSize = ceil(max(cell_landmarks) - min([cell_landmarks; 0, 0]));
    else
        imageSize = ceil(max(nuc_landmarks) - min([nuc_landmarks; 0, 0]));
    end
    % imageSize = flip(imageSize);
end

if isempty(cell_landmarks) && isempty(nuc_landmarks)
    error('cell landmarks and nuc landmarks must be provided at least one.');
end

if ~isempty(cell_landmarks)
    I_cell = poly2mask(cell_landmarks(:, 2), cell_landmarks(:, 1), imageSize(1), imageSize(2));
else
    I_cell = false(imageSize);
end

if ~isempty(nuc_landmarks)
    I_nuc = poly2mask(nuc_landmarks(:, 2), nuc_landmarks(:, 1), imageSize(1), imageSize(2));
else
    I_nuc = false(imageSize);
end
% convert landmarks to binary images. 
I = I_cell + I_nuc;

end