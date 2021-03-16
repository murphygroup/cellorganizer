function [hd, psnr, hd_mat, mbb_diagonal_len] = spharm_rpdm_reconstruction_error_computing(vertices_original, vertices_reconst)
% The function is calculate the reconstruction errors between the original
% and the reconstructed shape. 
% We calculate both the hausdorff distance as well as the peak
% signal-to-noise ratio. 
% hd_mat is an optional output for the directional hausdorff distances
% mbb_diagonal_len returns the length of the diagonal of the minimum
% bounding box

% Author: Xiongtao Ruan (xruan@andrew.cmu.edu)
%
% Copyright (C) 2013-2019 Murphy Lab
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


% calculate hausdorff distance
[hd, dhd, dhd_1] = hausdorff_distance_by_nearest_neighbor_search(vertices_original, vertices_reconst);
hd_mat = [dhd, dhd_1];

% calculate psnr
[rotmat, cornerpoints, volume, surface] = minboundbox(vertices_original(:, 1), vertices_original(:, 2), vertices_original(:, 3), 'v', 3);
mbb_diagonal_len = sqrt(sum((cornerpoints(1, :) - cornerpoints(7, :)) .^ 2));

psnr = 20 * log10(mbb_diagonal_len ./ hd);


end