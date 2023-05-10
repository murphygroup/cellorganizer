function [h] = spharm_rpdm_sample_or_reconstruct_images_figure(model, fileID, options)
% This function is used for illustration of shape evolution in the shape space for Spharm_RPDM method, 
% The shape evolution is performed through linear pathes in the shape space. 
% The two end points are chosen from the training image randomly. 

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
% 4/26/2023 R.F. Murphy Don't close open figures and don't make figure
% visible; add support for making movies (options.makemovie); ignore
% options.pair_method
% 5/1/2023 R.F. Murphy Switch video profile if deployed (Linux doesn't support MPEG-4)
% 5/2/2023 R.F. Murphy close figure window if making movie
% 5/3/2023 R.F. Murphy allow specification of movie profile through makemovie option
% 5/7/2023 R.F. Murphy pass fileID so that html can be written
% 5/11/2023 R.F. Murphy fix makemovie value comparisons

%close all;
f = figure('visible','off');

if nargin < 3
	options = struct();
end

rpdm_model = model.cellShapeModel;
train_score = rpdm_model.train_score;
train_coeff = rpdm_model.train_coeff;
mu = rpdm_model.mu;
scales = rpdm_model.scales;

options = ml_initparam(options, struct('num_steps', 10, 'makemovie', 'none'));

switch options.shape_evolution
   case 'random'
      % random sampling
      rand_inds = randperm(size(train_score, 1), 2);
      s_ind = rand_inds(1);
      t_ind = rand_inds(2);
   case 'first two shapes'
      s_ind = 1;
      t_ind = 2;
   case 'first and last shape'
      s_ind = 1;
      t_ind = size(train_score,1);
   case 'given'
      s_ind = options.pair_inds(1);
      t_ind = options.pair_inds(2);
      if any(options.pair_inds > size(train_score, 1)) || any(options.pair_inds <= 0)
         warning('The given pair %d and %d are not valid, using first two', options.pair_inds(1), options.pair_inds(2));
         s_ind = 1;
         t_ind = 2;
      end
   otherwise
      disp('Invalid shape_evolution option')
      return
end

s_yt = train_score(s_ind, :);
t_yt = train_score(t_ind, :);
steps = linspace(0, 1, options.num_steps + 1);
yt_evolve = s_yt .* (1 - steps') + t_yt .* steps';
sh_evolve = yt_evolve * train_coeff' + mu;

scale_evolve = scales(s_ind) .* (1 - steps) + scales(t_ind) .* steps;
scale_evolve = scale_evolve(:);

if any(strcmp(rpdm_model.components, 'nuc')) && any(strcmp(rpdm_model.components, 'cell'))
	sh_evolve = reshape(sh_evolve, size(sh_evolve, 1), [], 3, 2);
	sh_evolve = permute(sh_evolve, [2, 3, 1, 4]);
	cell_sh_evolve = sh_evolve(:, :, :, 1);
	nuc_sh_evolve = sh_evolve(:, :, :, 2);
    
	cell_center = rpdm_model.all_centers(1, :, :);
	nuc_center = rpdm_model.all_centers(2, :, :);
    switch rpdm_model.shape_model_type
        case 1        
			cell_centers_evolve = cell_center(:, :, s_ind) .* (1 - steps)' + cell_center(:, :, t_ind) .* steps';
			cell_centers_evolve = permute(cell_centers_evolve, [3, 2, 1]);
			nuc_centers_evolve = nuc_center(:, :, s_ind) .* (1 - steps)' + nuc_center(:, :, t_ind) .* steps';
			nuc_centers_evolve = permute(nuc_centers_evolve, [3, 2, 1]);

			cell_spharm_evolve = convert_3d_preshape_to_spharm(cell_sh_evolve, scale_evolve, cell_centers_evolve);
			nuc_spharm_evolve = convert_3d_preshape_to_spharm(nuc_sh_evolve, scale_evolve, nuc_centers_evolve);
		case 2
			cell_centers_evolve = cell_center(:, :, s_ind) .* (1 - steps)' + cell_center(:, :, t_ind) .* steps';
			cell_centers_evolve = permute(cell_centers_evolve, [3, 2, 1]);
			nuc_centers_evolve = nuc_center(:, :, s_ind) .* (1 - steps)' + nuc_center(:, :, t_ind) .* steps';
			nuc_centers_evolve = permute(nuc_centers_evolve, [3, 2, 1]);

			cell_spharm_evolve = convert_3d_preshape_to_spharm(cell_sh_evolve, scale_evolve, cell_centers_evolve);
			cell_spharm_evolve(1, :, :) = [];
			cell_spharm_evolve(1, :, :) = cell_spharm_evolve(1, :, :) + cell_centers_evolve;
			nuc_spharm_evolve = convert_3d_preshape_to_spharm(nuc_sh_evolve, scale_evolve, nuc_centers_evolve);
			nuc_spharm_evolve(1, :, :) = [];
			nuc_spharm_evolve(1, :, :) = nuc_spharm_evolve(1, :, :) + nuc_centers_evolve;
    end
else
    sh_evolve = reshape(sh_evolve, size(sh_evolve, 1), [], 3);
    sh_evolve = permute(sh_evolve, [2, 3, 1]);
	cur_center = rpdm_model.all_centers(1, :, :);
	centers_evolve = cur_center(:, :, s_ind) .* (1 - steps)' + cur_center(:, :, t_ind) .* steps';
    centers_evolve = permute(centers_evolve, [3, 2, 1]);    
	cur_spharm_evolve = convert_3d_preshape_to_spharm(sh_evolve, scale_evolve, centers_evolve);
	if any(strcmp(rpdm_model.components, 'nuc'))
		nuc_spharm_evolve = cur_spharm_evolve;
	elseif any(strcmp(rpdm_model.components, 'cell'))
		cell_spharm_evolve = cur_spharm_evolve;
	end
end

max_deg = rpdm_model.max_deg;
[vs, fs]=SpiralSampleSphere(4002);
% [vs fs] = sphereMesh([0 0 0 1]);
Zs = calculate_SPHARM_basis(vs, max_deg);


aspect_ration = 11/1.00;
if s_ind == 5 || t_ind == 5
	aspect_ration = 11/1;
end
figure_width = 1920;
figure_height = round(figure_width / aspect_ration);
dpi = 300;
set(gcf, 'Position', [-1920, 1, figure_width, figure_height]);
set(gcf, 'color', 'w')
%set(gcf, 'Visible', 'on')
border_width = 0.00;
border_hight = 0.05;
per_height = 1;
per_width = per_height / aspect_ration;

w_interval = -0.00;
h_interval = 0.00;

% define the bounding box for each figure, using cell shapes of source and
% targe cell. 
if any(strcmp(rpdm_model.components, 'cell')) 
    cell_vertices_s = real(Zs*cell_spharm_evolve(:, :, 1));
    cell_vertices_t = real(Zs*cell_spharm_evolve(:, :, end));
    cell_vertices_st = [cell_vertices_s; cell_vertices_t];
    bbox = [min(cell_vertices_st) - [2, 2, 0.5]; max(cell_vertices_st) + [2, 2, 0.5]];
elseif any(strcmp(rpdm_model.components, 'nuc')) 
    nuc_vertices_s = real(Zs*nuc_spharm_evolve(:, :, 1));
    nuc_vertices_t = real(Zs*nuc_spharm_evolve(:, :, end));
    nuc_vertices_st = [nuc_vertices_s; nuc_vertices_t];
    bbox = [min(nuc_vertices_st) - [2, 2, 0.5]; max(nuc_vertices_st) + [2, 2, 0.5]];
end

if ~strcmp(options.makemovie,'none')
    vidfile = VideoWriter('ShapeEvolutionMovie.avi',options.makemovie);
    open(vidfile);
end

for i = 1 : numel(steps)
	tic
	cell_vertices = [];
	nuc_vertices = [];	
	if any(strcmp(rpdm_model.components, 'cell')) 
		fvec_cell = cell_spharm_evolve(:, :, i);
		cell_vertices = real(Zs*fvec_cell);
		cell_faces = fs;
	end
	if any(strcmp(rpdm_model.components, 'nuc')) 
		fvec_nuc = nuc_spharm_evolve(:, :, i);
		Zvert_nuc = real(Zs*fvec_nuc);
	    [nuc_vertices, nuc_faces] = reorder_mesh_vertices(Zvert_nuc, fs, 'x-axis');
    end

    if strcmp(options.makemovie,'none')
    	curr_col = i;
        curr_row = 1;
        curr_position = [border_width + ( w_interval + per_width) * (curr_col - 1), 1 - (border_hight + (h_interval + per_height) * (curr_row -1) + per_height), per_width, per_height];
        axes('Position', curr_position);
    end
	
	if ~isempty(cell_vertices) && ~isempty(nuc_vertices)
		h = patch('vertices', cell_vertices, 'faces', cell_faces, 'FaceVertexCData', jet(size(cell_vertices,1)),'FaceColor','interp', 'EdgeColor', 'None', 'FaceAlpha', 'flat');
		alpha(0.5);
		hold on
		h = patch('vertices', nuc_vertices, 'faces', nuc_faces, 'FaceVertexCData', winter(size(nuc_vertices,1)),'FaceColor','interp', 'EdgeColor', 'None');
	elseif ~isempty(cell_vertices)
		h = patch('vertices', cell_vertices, 'faces', cell_faces, 'FaceVertexCData', jet(size(cell_vertices,1)),'FaceColor','interp', 'EdgeColor', 'None');
	elseif ~isempty(nuc_vertices)
		h = patch('vertices', nuc_vertices, 'faces', nuc_faces, 'FaceVertexCData', jet(size(nuc_vertices,1)),'FaceColor','interp', 'EdgeColor', 'None');
	end
	daspect([1 1 1])
    xlim(bbox(:, 1)')
    ylim(bbox(:, 2)')
    zlim(bbox(:, 3)')
	title(sprintf('%0.1f', steps(i)));
	axis off
	view([45, 45]);
	%camlight
	%lighting gouraud
	% axis tight
	set(findobj(gca, '-property', 'LineWidth'), 'LineWidth', 1.5);
	set(findobj(gca, '-property', 'fontsize'), 'fontsize', 15);
	set(findobj(gca, '-property', 'fontweight'), 'fontweight', 'bold');

    if ~strcmp(options.makemovie,'none')
        %saveas( gcf, ['frame' int2str(i) '.png'], 'png' );
        Frm=getframe(gcf);
        writeVideo(vidfile,Frm);
        hold off
    end

end

if ~strcmp(options.makemovie,'none')
    close(vidfile);
else
    saveas( f, 'show_shape_evolution.png', 'png' );
    fP = f.Position;
    f.Position = [fP(1) fP(2) round(fP(3)/2) round(fP(4)/2)];
    saveas( f, 'show_shape_evolution_thumbnail.png', 'png' ); %3/20/2023
    img2html(fileID, 'show_shape_evolution.png', ...
        'show_shape_evolution_thumbnail.png', ...
        ['Shape Evolution']);
end

close(f)

end



