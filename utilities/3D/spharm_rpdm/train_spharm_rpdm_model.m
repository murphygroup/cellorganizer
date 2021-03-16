function [spharm_rpdm_model] = train_spharm_rpdm_model(cell_params_fname, options)
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

% 09/19/2018 use spharm_rpdm to define the components, so that it is more
% flexible for only cell, only nuclear and joint shapes. 
% 02/07/2019 xruan also save the filenames of parameterization.
% 01/22/2021 R.F. Murphy save the final hausdorff distances into the model

options = ml_initparam(options, struct('latent_dim', 100, ...
                                       'maxDeg', 31));
                                   
spharm_rpdm = ml_initparam(options.spharm_rpdm, struct('shape_model_type', 1));
                                   
latent_dim = options.latent_dim;
max_deg = options.maxDeg; 
shape_model_type = spharm_rpdm.shape_model_type;

components = options.spharm_rpdm.components;
if numel(components) == 1
    rpdm_model_type = 'separate';
elseif any(strcmp(components, 'nuc')) && any(strcmp(components, 'cell')) 
    rpdm_model_type = 'joint';
end

nuc_spharm_descriptors = [];
if any(strcmp(components, 'nuc'))
    param_tmp = loadFiles(cell_params_fname, 'nuc');
    param_tmp = param_tmp(~cellfun(@isempty, param_tmp));
    % this line moves everything in param_tmp{i}.cell into param_temp{i}
    param_tmp = cellfun(@(x) x.nuc, param_tmp, 'Uniformoutput', false);
    nuc_spharm_descriptors = collect_descriptors_from_parameters(param_tmp, max_deg);
    nuc_spharm_distances = collect_distances_from_parameters(param_tmp);
end
cell_spharm_descriptors = [];
if any(strcmp(components, 'cell'))
    param_tmp = loadFiles(cell_params_fname, 'cell');
    param_tmp = param_tmp(~cellfun(@isempty, param_tmp));
    % this line moves everything in param_tmp{i}.cell into param_temp{i}
    param_tmp = cellfun(@(x) x.cell, param_tmp, 'Uniformoutput', false);    
    cell_spharm_descriptors = collect_descriptors_from_parameters(param_tmp, max_deg);
    cell_spharm_distances = collect_distances_from_parameters(param_tmp);
end

if strcmp(rpdm_model_type, 'separate')
    if exist(nuc_spharm_descriptors, 'var')
        all_spharm_descriptors = nuc_spharm_descriptors;
        all_spharm_distances = nuc_spharm_distances;
    else
        all_spharm_descriptors = cell_spharm_descriptors;
        all_spharm_distances = cell_spharm_distances;
    end
    all_centers = all_spharm_descriptors(1, :, :);
    all_spharm_descriptors_1 = all_spharm_descriptors;
    all_spharm_descriptors_1(1, :, :) = [];
    X_0 = [real(all_spharm_descriptors_1); imag(all_spharm_descriptors_1)];
    X = reshape(X_0, [], size(X_0, 3))';
elseif strcmp(rpdm_model_type, 'joint')
    all_spharm_descriptors = cat(1, cell_spharm_descriptors, nuc_spharm_descriptors);
    all_spharm_distances = cat(1, cell_spharm_distances, nuc_spharm_distances);
    switch shape_model_type
        case 1
            % remove center of the descriptors, cell and nuclear shapes are
            % centered by their own
            cell_centers = cell_spharm_descriptors(1, :, :);
            cell_spharm_descriptors(1, :, :) = [];
            nuc_centers = nuc_spharm_descriptors(1, :, :);
            nuc_spharm_descriptors(1, :, :) = [];

            X_train_cell_0 = [real(cell_spharm_descriptors); imag(cell_spharm_descriptors)];
            X_train_nuc_0 = [real(nuc_spharm_descriptors); imag(nuc_spharm_descriptors)];

            X_train_cell = reshape(X_train_cell_0, [], size(X_train_cell_0, 3))';
            X_train_nuc = reshape(X_train_nuc_0, [], size(X_train_nuc_0, 3))';
            X = [X_train_cell, X_train_nuc];            
        case 2
            % remove center of the descriptors, nuclear shapes are
            % centerized by the centers of corresponding cell shape
            cell_centers = cell_spharm_descriptors(1, :, :);
            cell_spharm_descriptors(1, :, :) = 0;
            nuc_spharm_descriptors(1, :, :) = nuc_spharm_descriptors(1, :, :) - cell_centers;
            nuc_centers = cell_centers;

            X_train_cell_0 = [real(cell_spharm_descriptors); imag(cell_spharm_descriptors)];
            X_train_nuc_0 = [real(nuc_spharm_descriptors); imag(nuc_spharm_descriptors)];

            X_train_cell = reshape(X_train_cell_0, [], size(X_train_cell_0, 3))';
            X_train_nuc = reshape(X_train_nuc_0, [], size(X_train_nuc_0, 3))';
            X = [X_train_cell, X_train_nuc];
    end
    all_centers = cat(1, cell_centers, nuc_centers);
end

[scales,coeff,score,latent,tsquared,explained,mu,train_score,train_explained,train_coeff] = CalculatePCA(X, latent_dim);

spharm_rpdm_model = struct();
spharm_rpdm_model.X = X;
spharm_rpdm_model.coeff = coeff;
spharm_rpdm_model.score = score;
spharm_rpdm_model.latent = latent;
spharm_rpdm_model.tsquared = tsquared;
spharm_rpdm_model.explained = explained;
spharm_rpdm_model.mu = mu;
spharm_rpdm_model.latent_dim = latent_dim;
spharm_rpdm_model.train_score = train_score;
spharm_rpdm_model.train_explained = train_explained;
spharm_rpdm_model.train_coeff = train_coeff;
spharm_rpdm_model.all_centers = all_centers;
spharm_rpdm_model.scales = scales;
spharm_rpdm_model.all_spharm_descriptors = all_spharm_descriptors;
spharm_rpdm_model.shape_model_type = shape_model_type;          % shape model type for different ways of joint/separate modeling
spharm_rpdm_model.rpdm_model_type = rpdm_model_type;            % joint or separate model
spharm_rpdm_model.components = components;                      % components
spharm_rpdm_model.max_deg = max_deg;
spharm_rpdm_model.hausdorff_distances = all_spharm_distances;

% other settings
spharm_rpdm_model.options = options;
spharm_rpdm_model.cell_params_fname = cell_params_fname;
spharm_rpdm_model.numimgs = size(cell_spharm_descriptors, 3);
spharm_rpdm_model.components = components;
spharm_rpdm_model.cell_params = cell_spharm_descriptors;
spharm_rpdm_model.nuc_params = nuc_spharm_descriptors;
spharm_rpdm_model.name = 'spharm-rpdm model of the cell';
spharm_rpdm_model.type = 'spharm_rpdm';
spharm_rpdm_model.version = 1.0;

end


function [all_spharm_descriptors] = collect_descriptors_from_parameters(input_params, max_deg)
% extract descriptors from the cell of structures of parameters. 

max_dim = (max_deg + 1)  ^ 2; 
num_params = numel(input_params);

all_spharm_descriptors = zeros((max_deg + 1) .^ 2, 3, num_params);

for i = 1 : num_params
    cur_fvec = input_params{i}.fvec;
    all_spharm_descriptors(1 : size(cur_fvec, 1), :, i) = cur_fvec;    
end

end

function [all_spharm_distances] = collect_distances_from_parameters(input_params)
% extract hausdorff distances from the cell of structures of parameters. 

num_params = numel(input_params);

for i = 1 : num_params
    all_spharm_distances(i) = input_params{i}.hd;    
end

end

