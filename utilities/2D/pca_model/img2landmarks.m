function [cur_landmarks] = img2landmarks(cur_image, pad_size, N_points)
% convert image to landmarks
% the parameterization method is (largely) based on the paper
% Pincus, Zachary, and J. A. Theriot. "Comparison of quantitative methods for cell‐shape analysis." Journal of microscopy 227.2 (2007): 140-156.

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


if nargin < 2
	pad_size = 3;
end
if nargin < 3
	N_points = 1000;
end

if ~any(cur_image(:) > 0) 
    warning('The image is empty, just return a NaN landmark')
    cur_landmarks = NaN(N_points, 2);
    return;
end
CC = bwconncomp(cur_image > 0);
if CC.NumObjects > 1
    pixel_nums = cellfun(@numel, CC.PixelIdxList);
    [~, max_ind] = max(pixel_nums);
    max_PixelIdxs = CC.PixelIdxList{max_ind};
    BW_1 = false(size(cur_image));
    BW_1(max_PixelIdxs) = true;
    cur_image = BW_1;
end
cur_cell_boundary = bwboundaries(cur_image > 0);  
cur_cell_boundary = cur_cell_boundary{1};

% direct interpolation
cur_cell_boundary = padarray(cur_cell_boundary(1 : end - 1, :), [pad_size, 0], 'circular', 'post');

interp_interval = 0.01;
interp_i = 1 : interp_interval : size(cur_cell_boundary, 1);
interp_y = interp1(cur_cell_boundary(:, 1), interp_i, 'spline');
interp_x = interp1(cur_cell_boundary(:, 2), interp_i, 'spline');
if false
    figure, scatter(interp_x, interp_y);
end
ind = find((interp_y == cur_cell_boundary(1, 1)) .* (interp_x == cur_cell_boundary(1, 2)) == 1, 1, 'last');
interp_BW = [interp_y(1 : ind)', interp_x(1 : ind)'];

% compute circumstance
cur_intervals_lengths = sqrt(sum(diff(interp_BW) .^ 2, 2));
cur_circumstance = sum(cur_intervals_lengths);
per_length = cur_circumstance / N_points;


% find sampling points. Just start from the first point for simple. 
cur_sampling_points = zeros(N_points, 2);
cur_sampling_points(1, :) = interp_BW(1, :);
accumulate_dist = 0;
N_add = 1;
for j = 2 : size(interp_BW, 1)
    accumlate_dist_old = accumulate_dist;
    accumulate_dist = accumulate_dist + cur_intervals_lengths(j - 1);
    if accumulate_dist / per_length >= N_add
        if accumlate_dist_old / per_length < N_add
            d1 = N_add * per_length - accumlate_dist_old;
            d2 = accumulate_dist - N_add * per_length;
            cur_sampling_points(N_add + 1, :) = (d2 * interp_BW(j - 1, :) + d1 * interp_BW(j, :)) / (d1 + d2);
            N_add = N_add + 1;                
        else
            continue;
        end
    end

end
if false
    close all;
    figure, scatter(interp_BW(:, 1), interp_BW(:, 2))
    figure, scatter(cur_sampling_points(:, 1), cur_sampling_points(:, 2));
end
% centering using the landmarks
cur_sampling_points = cur_sampling_points(1 : N_points, :);
% 01/01/2018
% cur_landmarks = cur_sampling_points - mean(cur_sampling_points);
cur_landmarks = cur_sampling_points; 

end
