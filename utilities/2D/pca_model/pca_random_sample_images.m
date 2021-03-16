function [nucimg, cellimg] = pca_random_sample_images( model, options )
% randomly sample images, idea: first randomly pick two point in the shape space, and then 
% randomly sample a point uniformly lay in the line segment. At last reconstruct the landmarmks
% of this point. 

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



pca_model = model.cellShapeModel;
train_score = pca_model.train_score;
train_coeff = pca_model.train_coeff;
mu = pca_model.mu;

% random sampling
rand_inds = randperm(size(train_score, 1), 2);
train_point_1 = train_score(rand_inds(1), :);
train_point_2 = train_score(rand_inds(2), :);
lambda = rand();
rand_point = train_point_1 * lambda + train_point_2 * (1 - lambda);

reconst_landmarks = rand_point * train_coeff' + mu;
reconst_landmarks = reshape(reconst_landmarks', [], 2); 
centers = pca_model.centers;
scales = pca_model.scales;
center_1 = centers(:, :, rand_inds(1));        
center_2 = centers(:, :, rand_inds(2)); 
cur_center = center_1 * lambda +center_2 * (1 - lambda);

switch pca_model.options.pca_options.shape_type
    case {'shape', 'preshape'}
        scale_1 = scales(rand_inds(1));
        scale_2 = scales(rand_inds(2));
        cur_scale = scale_1 * lambda +scale_2 * (1 - lambda);
        reconst_landmarks = reconst_landmarks .* cur_scale + cur_center;
    case {'shape_size'}
        reconst_landmarks = reconst_landmarks + cur_center;        
end

cell_reconst_landmarks = [];
nuc_reconst_landmarks = [];
if any(strcmp(pca_model.components, 'nuc')) && any(strcmp(pca_model.components, 'cell')) 
    cell_reconst_landmarks = reconst_landmarks(1 : 1000, :);
    nuc_reconst_landmarks = reconst_landmarks(1001 : 2000, :);
elseif any(strcmp(pca_model.components, 'nuc'))
    nuc_reconst_landmarks = reconst_landmarks(1 : 1000, :);
elseif any(strcmp(pca_model.components, 'cell'))
    cell_reconst_landmarks = reconst_landmarks(1 : 1000, :);
end

% convert landmarks to binary images. 
imageSize = options.imageSize;
cellimg = [];
nucimg = [];

if ~isempty(cell_reconst_landmarks)
    cellimg = poly2mask(cell_reconst_landmarks(:, 2), cell_reconst_landmarks(:, 1), imageSize(1), imageSize(2));
end

if ~isempty(nuc_reconst_landmarks)
    nucimg = poly2mask(nuc_reconst_landmarks(:, 2), nuc_reconst_landmarks(:, 1), imageSize(1), imageSize(2));
end

cellimg = convert_binary_image_to_outline(cellimg);
nucimg = convert_binary_image_to_outline(nucimg);



end



