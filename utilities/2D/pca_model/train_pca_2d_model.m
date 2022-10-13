function [pca_model] = train_pca_2d_model(cell_params_fname, options)
% take in binary mask cell shape, use some smoothing methods
% (interpolation, spline smoothing, kernel smoothing) to smooth the cell
% boundary, and then resample the cell boundary points, evening. After
% that do translation to centering the cell shapes. 

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

% 03/23/2018 xruan minor fix of saving model content
% 01/13/2019 xruan add option of normalization (scale, orientation), and set
% the normalization of orientation as default options (use shape and scale 
% for single component); for joint modeling, use cell alignment with size as default). 



options = ml_initparam(options, struct('latent_dim', 10));
latent_dim = options.latent_dim;

if ~isfield(options, 'pca_options')
    options.pca_options = struct();
end

% set default alignment method for joint shape, and the input shape type
% for pca model. The alignment method only valid for shape size and shape
pca_options = ml_initparam(options.pca_options, struct('shape_type', 'shape_size', ...  % other choices: shape, preshape, original
                                                       'alignment_method', 'cell_align', ...  % other choices: separate_align, joint align
                                                       'outline_scaling_method', 'joint_scaling')); % other choices: separate_scaling

components = options.components;
if numel(components) == 1
    pca_model_type = 'separate';
elseif any(strcmp(components, 'nuc')) && any(strcmp(components, 'cell')) 
    pca_model_type = 'joint';
end

nuc_params = [];
if any(strcmp(components, 'nuc'))
    param_tmp = loadFiles(cell_params_fname, 'nuc');
    param_tmp = param_tmp(~cellfun(@isempty, param_tmp));
    nuc_params = cellfun(@(x) x.nuc, param_tmp, 'UniformOutput', false);
    nuc_landmarks = cat(3, nuc_params{:});
end
cell_params = [];
if any(strcmp(components, 'cell'))
    param_tmp = loadFiles(cell_params_fname, 'cell');
    param_tmp = param_tmp(~cellfun(@isempty, param_tmp));
    cell_params = cellfun(@(x) x.cell, param_tmp, 'UniformOutput', false);
    cell_landmarks = cat(3, cell_params{:});
end

if strcmp(pca_model_type, 'separate')
    if exist('cell_landmarks', 'var')
        all_landmarks = cell_landmarks;
    else
        all_landmarks = nuc_landmarks;
    end
    % postprocess
    [all_landmarks_shape, scales, centers, angles] = pca_shape_processing({all_landmarks}, pca_options);
elseif strcmp(pca_model_type, 'joint')
    all_landmarks = cat(1, cell_landmarks, nuc_landmarks); 
    % postprocess
    all_landmarks_cell = {cell_landmarks, nuc_landmarks};
    [all_landmarks_shape, scales, centers, angles, cell_nuc_scales, cell_nuc_centers, cell_nuc_angles] = pca_shape_processing(all_landmarks_cell, pca_options);
end

X = reshape(all_landmarks_shape, [], size(all_landmarks_shape, 3))';
switch pca_options.shape_type
    case {'original', 'preshape', 'shape'}
    case 'shape_size'
        X =  X .* squeeze(scales);
end
        
% X = X - mean(X);

[coeff,score,latent,tsquared,explained, mu] = pca(X, 'NumComponents', latent_dim);

disp(size(score))

train_score = score;
train_explained = sum(explained);
train_coeff = coeff;


pca_model = struct();
pca_model.X = X;
pca_model.coeff = coeff;
pca_model.score = score;
pca_model.latent = latent;
pca_model.tsquared = tsquared;
pca_model.explained = explained;
pca_model.mu = mu;
pca_model.latent_dim = latent_dim;
pca_model.train_score = train_score;
pca_model.train_explained = train_explained;
pca_model.train_coeff = train_coeff;
pca_model.all_landmarks = all_landmarks;
pca_model.all_landmarks_shape = all_landmarks_shape;
pca_model.scales = scales;
pca_model.centers = centers;
pca_model.angles = angles;
if strcmp(pca_model_type, 'joint')
    pca_model.cell_nuc_scales = cell_nuc_scales;
    pca_model.cell_nuc_centers = cell_nuc_centers;
    pca_model.cell_nuc_angles = cell_nuc_angles;
end

% other settings
options.pca_options = pca_options;
pca_model.options = options;
if ~isempty(cell_params)
    pca_model.numimgs = numel(cell_params);
else
    pca_model.numimgs = numel(nuc_params);
end
pca_model.components = components;
pca_model.cell_params = cell_params;
pca_model.nuc_params = nuc_params;

% if use shape or shape_size the cell and nuclear params are aligne outlines
switch pca_options.shape_type
    case {'original', 'preshape'}
        pca_model.cell_params = cell_params;
        pca_model.nuc_params = nuc_params;
    case {'shape', 'shape_size'}
        pca_model.orig_cell_params = cell_params;
        pca_model.orig_nuc_params = nuc_params;
        all_params = all_landmarks_shape .* scales + centers;
        
        if strcmp(pca_model_type, 'separate')
            if any(strcmp(components, 'cell'))
                pca_model.cell_params = mat2cell(all_params, size(all_params, 1), size(all_params, 2), ones(size(all_params, 3), 1));
            else
                pca_model.nuc_params = mat2cell(all_params, size(all_params, 1), size(all_params, 2), ones(size(all_params, 3), 1));
            end
        elseif strcmp(pca_model_type, 'joint')
            all_cell_params = all_params(1 : size(cell_params{1}, 1), :, :);
            all_nuc_params = all_params(size(cell_params{1}, 1) + 1 : end, :, :);
            pca_model.cell_params = mat2cell(all_cell_params, size(all_cell_params, 1), size(all_cell_params, 2), ones(size(all_cell_params, 3), 1));
            pca_model.nuc_params = mat2cell(all_nuc_params, size(all_nuc_params, 1), size(all_nuc_params, 2), ones(size(all_nuc_params, 3), 1));
        end       
end
        
pca_model.name = 'pca model of the cell';
pca_model.type = 'pca';
pca_model.version = 1.0;

end


