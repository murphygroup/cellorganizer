function [param] = show_SPHARM_RPDM_Shape_Space_Figure(model, labels, param)
% SHOWPCASHAPESPACEFIGURE 

% Xiongtao Ruan
%
% Copyright (C) 2015-2018 Murphy Lab
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
% For additional information visit http://murphylab.web.cmu.edu/ or
% send email to murphy@cmu.edu

% 03/23/2018 copied from showPCAShapeSpaceFigure.m
% 05/24/2018 add support for labeling
% 02/26/2019 make enhancement to throw errors when the model is invalid. 
% 03/05/2021 R. F. Murphy fix hausdorff cutoff to consider both nuc and cell
% 03/23/2023 R. F. Murphy fix latent dimension check

if ~exist('param', 'var')
    param = [];
end

train_score = model.train_score;
nimgs = size(train_score, 1);

if size(train_score, 2) > 3
    [~, synthorder] = sort(train_score(:,4), 'Ascend');
else
    synthorder = 1:nimgs;
end

param = ml_initparam(param, struct( ...
    'skipmissing', true, ...
    'rebuild_pos', false, ...
    'plot_dims', [1,2,3], ...
    'cm', @jet, ...
    'traces', [], ...
    'embedfctn', @cmdscale, ...
    'subsize', 400, ...
    'figure_scale_ratio', 0.5, ...
    'synthorder', synthorder, ...
    'rotate', 0 ...
    ));
cellnums = 1:nimgs;

if ~exist('labels', 'var') labels = []; end

if strcmpi(labels,'unique')
    labels = [1:nimgs];
elseif isempty(labels)
    labels = ones(nimgs, 1);
end

skipmissing = param.skipmissing;
cm = param.cm;
plot_dims = param.plot_dims;
rebuild_pos = param.rebuild_pos;

% check whether the model is valid
if size(train_score, 1) < 5 %3/23/2023
    error('The minimum number of cells in the model is 5, while the current model only contains %d cell(s). Please increase the number of cells', size(train_score, 1));
end

latent_dims = size(train_score, 2);
if latent_dims < max(plot_dims)
    error('At least one of chosen plot dimensions [ %s ] is larger than the latent dimension of the model %d. Please reset param.plot_dims ', num2str(plot_dims), latent_dims);
end


% keepinds = ~any(isnan(train_score),2);
if isfield(model, 'hausdorff_distances') && isfield(param, 'hd_threshold')
    keepinds = ~(any(isnan(train_score),2) | (model.hausdorff_distances>param.hd_threshold)');
    %     labels = labels(keepinds);
    if strcmp(model.rpdm_model_type,'joint')
        keepinds = keepinds(:,1) & keepinds(:,2);
    end    
    cellnums = cellnums(keepinds);
    %     embed_pos = model.positions(keepinds,:);
end

[ulabels, ~, labelinds ] = unique(labels);

labels = labelinds;

if strcmpi(class(cm), 'function_handle') && length(ulabels) > 1
    colors = cm(length(ulabels))*0.8;
else
    colors = cm;
end

param.colors = colors;
set(gcf, 'color', 'w')
hold on


normalize_coords = train_score;
normalize_coords = normalize_coords ./ max(abs(normalize_coords)) * param.subsize * param.figure_scale_ratio;
if numel(plot_dims) < 3
    normalize_coords = [normalize_coords(:, plot_dims), zeros(size(normalize_coords, 1), 1)];
else
    normalize_coords = normalize_coords(:, plot_dims);
end

% construct mesh from descriptors
all_spharm_descriptors = model.all_spharm_descriptors;
if numel(model.components) > 1
    all_spharm_descriptors = reshape(all_spharm_descriptors, [], 2, size(all_spharm_descriptors, 2), size(all_spharm_descriptors, 3));
    all_spharm_descriptors = permute(all_spharm_descriptors, [1, 3, 2, 4]);
    mesh_func = @(x) reconstruct_spharm_descriptor_to_mesh(all_spharm_descriptors(:, :, :, x), model.components);
else
    mesh_func = @(x) reconstruct_spharm_descriptor_to_mesh(all_spharm_descriptors(:, :, x), model.components);
end

param.synthorder = param.synthorder(ismember(param.synthorder, cellnums));
for imnum = 1:length(param.synthorder)
    
    cellind = param.synthorder(imnum);
%     try
        cellimg = mesh_func(cellind);
        
        if isempty(cellimg) %|| any(cellimg(:) .* cellbounds(:))
            continue;
        end
%        if strcmpi(class(cm), 'function_handle')
        if length(ulabels) > 1
            param.colormap = param.colors(labels(cellind),:);
        else
            param.colormap = param.colors;
        end
        %center the image on the point where it belongs
        place_cell_nuc_3d_mesh(cellimg, normalize_coords(cellind, :), param);
        hold on
%     catch err
%         disp(['Unable to plot image: ' num2str(imnum)] );
%         getReport( err )
%     end
end

axis equal
axis tight

param.axis = axis;
hold off
disp(['Number of objects: ' num2str(length(param.synthorder))]);
%
% embed_final = nan(nimgs, size(embed_pos,2));
% embed_final(cellnums,:) = embed_pos;
%
% locs_final = nan(nimgs, 2);
% locs_final(cellnums,:) = locs;
%
% param.cellnums = cellnums;

end


function [colormap_func] = get_colormap(ind)

colormaps_cell = {'jet', 'parula', 'hsv', 'hot', 'cool', 'spring', 'summer', 'autumn', 'winter', 'gray', 'bone', 'copper', 'pink'};

ind = rem(ind, numel(colormaps_cell));
if ind == 0 
    ind = numel(colormaps_cell);
end

colormap_func = str2func(colormaps_cell{ind});

end
