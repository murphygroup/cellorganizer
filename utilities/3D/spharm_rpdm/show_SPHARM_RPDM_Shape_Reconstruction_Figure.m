function [figure_handle_cell, reconstruction_error_mat, image_inds] = show_SPHARM_RPDM_Shape_Reconstruction_Figure(model, options )
% This function is used for illustration of shape reconstuction for
% Spharm_RPDM method for given training image. Here, we visualize the
% figure, as well as reporting the reconstruction errors: hausdorff
% distance and peak signal-to-noise ratio
% The first output is for the figure handles
% The second output is to output the reconstruction errors as a mat file,
% which is optional. 
% The third output is the image indices in the plot

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

close all;

if nargin < 2
	options = struct();
end

options = ml_initparam(options, struct( ...
    'do_rotate_back', true, ... % decide whether to rotate back the images if there is alignment
    'image_inds', [1:9], ...              % image_inds is a list of integers indicating the train image to be reconstructed. 
    'parameterization_location', '', ...    % location of parameterization, by default it is provide in the model, in case not, manually provide it here.
    'figure_location', pwd,  ... % The location to save the figures
    'dpi', 300 ...  % resolution for the figure
    ));

rpdm_model = model.cellShapeModel;
train_score = rpdm_model.train_score;
train_coeff = rpdm_model.train_coeff;
mu = rpdm_model.mu;
scales = rpdm_model.scales;

do_rotate_back = options.do_rotate_back;
image_inds = options.image_inds;
parameterization_location = options.parameterization_location;
figure_location = options.figure_location;
dpi = options.dpi;

% check whether the parameterization files exists
cell_params_fname = rpdm_model.cell_params_fname;

cell_params_fname_to_plot = cell_params_fname(image_inds);
for i = 1 : numel(cell_params_fname_to_plot)
    cur_cell_param_fname = cell_params_fname_to_plot{i};
    if ~exist(cur_cell_param_fname, 'file')
        if isempty(parameterization_location)
            error('The original parameterization for image %d does not exist! \nPlease provide the directory location for the parameterizations (options.parameterization_location)', image_inds(i)); 
        else
            cell_params_fname_to_plot{i} = sprintf('%s/param%d.mat', image_inds(i));
        end        
    end
end

max_deg = rpdm_model.max_deg;
[vs, fs]=SpiralSampleSphere(4002);
% [vs fs] = sphereMesh([0 0 0 1]);
Zs = calculate_SPHARM_basis(vs, max_deg);

reconstruction_error_mat = zeros(numel(image_inds), 2);
figure_handle_cell = cell(numel(image_inds), 2);

for i = 1 : numel(image_inds)
    cur_image_ind = image_inds(i);
    cur_cell_params_fname = cell_params_fname_to_plot{i};
    a = load(cur_cell_params_fname);
    
    if any(strcmp(rpdm_model.components, 'cell')) 
        vertices_cell = a.cell.vertices;
        faces_cell = a.cell.faces;
        if do_rotate_back && a.cell.postprocess
            cell_R = a.cell.R;
            cell_rotation_center = a.cell.rotation_center;
            vertices_cell = (vertices_cell - cell_rotation_center) * cell_R + cell_rotation_center;            
        end        
    end
    
    if any(strcmp(rpdm_model.components, 'nuc'))
        vertices_nuc = a.nuc.vertices;
        faces_nuc = a.nuc.faces;
        if do_rotate_back && a.nuc.postprocess
            nuc_R = a.nuc.R;
            nuc_rotation_center = a.nuc.rotation_center;
            vertices_nuc = (vertices_nuc - nuc_rotation_center) * nuc_R + nuc_rotation_center;
        end
    end
    
    % reconstruct 
    train_score = rpdm_model.train_score;
    cur_train_score = train_score(cur_image_ind, :);
    cur_center = rpdm_model.all_centers(:, :, cur_image_ind);
    cur_scale = rpdm_model.scales(cur_image_ind, :);
    reconst_spharm_descriptors = cur_train_score * train_coeff' + mu;

    cell_spharm_descriptor = [];
    nuc_spharm_descriptor = [];
    if any(strcmp(rpdm_model.components, 'nuc')) && any(strcmp(rpdm_model.components, 'cell')) 
        reconst_spharm_descriptors = reshape(reconst_spharm_descriptors', [], 3, 2);     
        cell_reconst_descriptors = reconst_spharm_descriptors(:, :, 1);
        nuc_reconst_descriptors = reconst_spharm_descriptors(:, :, 2);
        cell_center = cur_center(1, :);
        nuc_center = cur_center(2, :);
        switch rpdm_model.shape_model_type
            case 1 
                cell_spharm_descriptor = convert_3d_preshape_to_spharm(cell_reconst_descriptors, cur_scale, cell_center);
                nuc_spharm_descriptor = convert_3d_preshape_to_spharm(nuc_reconst_descriptors, cur_scale, nuc_center);
            case 2
                cell_spharm_descriptor = convert_3d_preshape_to_spharm(cell_reconst_descriptors, cur_scale, cell_center);
                cell_spharm_descriptor(1, :) = [];
                cell_spharm_descriptor(1, :) = cell_spharm_descriptor(1, :) + cell_center;
                nuc_spharm_descriptor = convert_3d_preshape_to_spharm(nuc_reconst_descriptors, cur_scale, nuc_center);
                nuc_spharm_descriptor(1, :) = [];
                nuc_spharm_descriptor(1, :) = nuc_spharm_descriptor(1, :) + nuc_center;            
        end
    elseif any(strcmp(rpdm_model.components, 'nuc'))
        nuc_reconst_descriptors = reconst_spharm_descriptors;
        nuc_spharm_descriptor = convert_3d_preshape_to_spharm(nuc_reconst_descriptors, cur_scale, cur_center);
    elseif any(strcmp(rpdm_model.components, 'cell'))
        cell_reconst_descriptors = reconst_spharm_descriptors;
        cell_spharm_descriptor = convert_3d_preshape_to_spharm(cell_reconst_descriptors, cur_scale, cur_center);    
    end
    
    Zvert_cell = [];
    if ~isempty(cell_spharm_descriptor)
        fvec_cell = cell_spharm_descriptor;
        Zvert_cell = real(Zs * fvec_cell);
        if do_rotate_back
            Zvert_cell = (Zvert_cell - cell_rotation_center) * cell_R + cell_rotation_center;                        
        end
    end
    
    Zvert_nuc = [];
    if ~isempty(nuc_spharm_descriptor)
        fvec_nuc = nuc_spharm_descriptor;
        Zvert_nuc = real(Zs * fvec_nuc);
        if do_rotate_back
            Zvert_nuc = (Zvert_nuc - nuc_rotation_center) * nuc_R + nuc_rotation_center;            
        end        
    end
    
    % calculate reconstruction errors
    cell_reconstruction_error = [];
    nuc_reconstruction_error = [];
    all_reconstruction_error = [];
    if ~isempty(Zvert_cell)
        [cell_hd, cell_psnr] = spharm_rpdm_reconstruction_error_computing(vertices_cell, Zvert_cell);
        cell_reconstruction_error = [cell_hd, cell_psnr];
        all_reconstruction_error = cell_reconstruction_error;        
    end
    
    if ~isempty(Zvert_nuc)
        [nuc_hd, nuc_psnr] = spharm_rpdm_reconstruction_error_computing(vertices_nuc, Zvert_nuc);
        nuc_reconstruction_error = [nuc_hd, nuc_psnr];
        all_reconstruction_error = nuc_reconstruction_error;        
    end
    
    % if joint modeling, use average error of cell and nuclear errors
    if ~isempty(Zvert_cell) && ~isempty(Zvert_nuc)
        all_reconstruction_error = (cell_reconstruction_error + nuc_reconstruction_error) / 2;
    end
    reconstruction_error_mat(i, :) = all_reconstruction_error;
    
    % plot the reconstructions
    figure, 
    aspect_ration = 10/6;
    figure_width = 1920;
    figure_height = round(figure_width / aspect_ration);
    border_width = 0;
    border_hight = 0.00;
    per_height = 0.9;
    per_width = per_height / aspect_ration;
    w_interval = -0.05;
    h_interval = -0.1;

    set(gcf, 'Position', [0, 0, figure_width, figure_height]);
    set(gcf, 'color', 'w')
    set(gcf, 'Visible', 'on')
    axis off;
    
    title_string = {'Original', 'Reconstruction'};
    for j = 1 : 2
        if j == 1
            if ~isempty(vertices_cell)
                vertices_cell_j = vertices_cell;
                faces_cell_j = faces_cell;
                [vertices_cell_j, faces_cell_j] = reorder_mesh_vertices(vertices_cell_j, faces_cell_j, 'given-direction', [1, 0, 0]);                
            end
            if ~isempty(vertices_nuc)
                vertices_nuc_j = vertices_nuc;
                faces_nuc_j = faces_nuc;
                [vertices_nuc_j, faces_nuc_j] = reorder_mesh_vertices(vertices_nuc_j, faces_nuc_j, 'given-direction', [1, 0, 0]);
            end            
        elseif j == 2
            if ~isempty(Zvert_cell)
                vertices_cell_j = Zvert_cell;
                faces_cell_j = fs;
                [vertices_cell_j, faces_cell_j] = reorder_mesh_vertices(vertices_cell_j, faces_cell_j, 'given-direction', [1, 0, 0]);                
            end
            if ~isempty(Zvert_nuc)
                vertices_nuc_j = Zvert_nuc;
                faces_nuc_j = fs;
                [vertices_nuc_j, faces_nuc_j] = reorder_mesh_vertices(vertices_nuc_j, faces_nuc_j, 'given-direction', [1, 0, 0]);
            end            
        end
        curr_col = j;
        curr_row = 1;
        curr_position = [border_width + ( w_interval + per_width) * (curr_col - 1), 1 - (border_hight + (h_interval + per_height) * (curr_row -1) + per_height), per_width, per_height];
        h = axes('Position', curr_position);
        figure_handle_cell{i, j} = h;
        
        if ~isempty(vertices_cell_j) && ~isempty(vertices_nuc_j)
            patch('vertices', vertices_cell_j, 'faces', faces_cell_j, 'FaceVertexCData', jet(size(vertices_cell_j,1)),'FaceColor','interp', 'EdgeColor', 'None', 'FaceAlpha', 'flat');
            alpha(0.5);
            hold on
            patch('vertices', vertices_nuc_j, 'faces', faces_nuc_j, 'FaceVertexCData', winter(size(vertices_nuc_j,1)),'FaceColor','interp', 'EdgeColor', 'None');
        elseif ~isempty(vertices_cell_j)
            patch('vertices', vertices_cell_j, 'faces', faces_cell_j, 'FaceVertexCData', jet(size(vertices_cell_j,1)),'FaceColor','interp', 'EdgeColor', 'None');
        elseif ~isempty(vertices_nuc_j)
            patch('vertices', vertices_nuc_j, 'faces', faces_nuc_j, 'FaceVertexCData', jet(size(vertices_nuc_j,1)),'FaceColor','interp', 'EdgeColor', 'None');
        end
        daspect([1 1 1])
        title(title_string{j}, 'Units', 'normalized', 'position', [0.5, 0.9, 0], 'FontSize', 21);
        axis off
        view([45, 45]);
        camlight
        lighting gouraud
        if j == 2
            ax = gca;
            label_string = {sprintf('HD: %0.2f', all_reconstruction_error(1)), sprintf('PSNR: %0.2f dB', all_reconstruction_error(2))};            
            xlabel(ax, label_string, 'HorizontalAlignment','center', 'Units', 'normalized', 'Position', [0.5, 0.10, 0], 'FontSize', 14)
            ax.XLabel.Visible = 'on';
        end
    end
    set(findobj(gca, '-property', 'fontweight'), 'fontweight', 'bold');
    figure_filename = sprintf('%s/image_%d_reconstruction_illustration.png', figure_location, cur_image_ind);
    export_fig(figure_filename, '-opengl', '-png', '-a1', '-nocrop', ['-r', num2str(dpi)]);
end


end

