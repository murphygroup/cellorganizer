function result = ml_gencellshape3d( model, nucleus, param )
% ML_GENCELLSHAPE3D generates a cell shape from the eigen shape model of the
% radius ratio of the nuclei to the cell membrane.

% Author: Ivan E. Cao-Berg
%
% Copyright (C) 2012-2019 Murphy Lab
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

% Created April 7, 2012 R.F. Murphy from tp_gencellshape by T. Peng
%
% April 7, 2012 R.F. Murphy Use param.samp_rate consistently as number
% of angles in model
% April 8, 2012 R.F. Murphy Add models for cell/nuclear height ratio
% and position of nucleus
% April 9, 2012 R.F. Murphy Regenerate and fill nuclear surface to resolve
% registration issues
% July 25, 2012 R.F. Murphy Handle case where nuclear bottom position
% distribution has zero standard deviation
% June 29, 2013 R.F. Murphy Add param.generatemeanshape to get mean shape
% July 1, 2013 R.F. Murphy Add default param.generatemeanshape=false

if nargin < 2
    error('2 or 3 arguments are required');
end

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param,...
    struct('xsize',1024,...
    'ysize',1024,...
    'samp_rate',360,...
    'alpha',0.2,...
    'generatemeanshape',false,...
    'synthesis','all'));

nucimgsize = nucleus.nucimgsize;

% should we use actual center of mass of nucleus here?
xcenter = nucimgsize(1) / 2;
ycenter = nucimgsize(2) / 2;

norm_std.name = 'mvn';
norm_std.mu = zeros(1,size(model.modeShape.const,2));
norm_std.sigma = eye(size(model.modeShape.const,2));

if param.generatemeanshape
    shape_vec = model.meanShape.const;
else
    syn_score = ml_rnd(norm_std);
    while tp_pvalue(syn_score,norm_std) < param.alpha
        syn_score = ml_rnd(norm_std);
    end
    syn_score = syn_score .* sqrt(model.eigens.stat');

    shape_vec = model.meanShape.const + (model.modeShape.const * syn_score')';
end

instanceheight = length(shape_vec)/param.samp_rate;
nuc2cell_ratio = reshape(shape_vec,[instanceheight param.samp_rate]);

% Up-sampling to get correct height of cell
% model contains Gaussian model for ratio of cell height to nuclear height
% draw instance from it and use it to find final cell height
try
    f_cellnucratio = model.cellnucratio.stat;
    cell_nuc_ratio = ml_rnd(f_cellnucratio);
    while tp_pvalue(cell_nuc_ratio,f_cellnucratio) < param.alpha
        cell_nuc_ratio = ml_rnd(f_cellnucratio);
    end
catch
    cell_nuc_ratio = (nucimgsize(3)+2)/nucimgsize(3);
end

cell_height = round(nucimgsize(3)*cell_nuc_ratio);
if cell_height~=instanceheight
    t = 0:(instanceheight-1)/cell_height:instanceheight-1;
    nuc2cell_ratio_interp = zeros(cell_height+1,param.samp_rate);
    for i = 1:param.samp_rate
        nuc2cell_ratio_interp(:,i) = interp1(0:instanceheight-1,nuc2cell_ratio(:,i),t);
    end
else % just copy if already the correct height
    nuc2cell_ratio_interp = nuc2cell_ratio;
end
nuc2cell_ratio_interp(:,end+1) = nuc2cell_ratio_interp(:,1);

delta = 2*pi/param.samp_rate;
Phi = -pi:delta:pi;

cellimg = zeros(nucimgsize(1),nucimgsize(2),cell_height);
nucimg = cellimg;

% Create mesh
n_angles = length(Phi)-1;
n_slices = cell_height;
cellmesh = struct('vertices', zeros(n_slices * n_angles + 2, 3), 'faces', zeros((n_slices-1) * n_angles * 2 + n_angles * 2, 3));
% Bottom
bottom_face_indices = (1:n_angles) + (n_slices-1) * n_angles * 2;
bottom_center_vertex_index = size(cellmesh.vertices, 1) - 1;
cellmesh.faces(bottom_face_indices, 1) = 1:n_angles;
cellmesh.faces(bottom_face_indices, 2) = circshift(cellmesh.faces(bottom_face_indices, 1), 1);
cellmesh.faces(bottom_face_indices, 3) = bottom_center_vertex_index;
cellmesh.vertices(bottom_center_vertex_index, 1) = xcenter;
cellmesh.vertices(bottom_center_vertex_index, 2) = ycenter;
cellmesh.vertices(bottom_center_vertex_index, 3) = 1;
% Top
top_face_indices = (1:n_angles) + (n_slices-1) * n_angles * 2 + n_angles;
top_center_vertex_index = size(cellmesh.vertices, 1);
cellmesh.faces(top_face_indices, 1) = (1:n_angles) + (n_slices-1) * n_angles;
cellmesh.faces(top_face_indices, 2) = circshift(cellmesh.faces(top_face_indices, 1), -1);
cellmesh.faces(top_face_indices, 3) = top_center_vertex_index;
cellmesh.vertices(top_center_vertex_index, 1) = xcenter;
cellmesh.vertices(top_center_vertex_index, 2) = ycenter;
cellmesh.vertices(top_center_vertex_index, 3) = n_slices;
% Sides
for i = 1:n_slices-1
    side_face_indices1 = (1:n_angles) + (i-1) * n_angles * 2;
    cellmesh.faces(side_face_indices1, 1) = (1:n_angles) + (i-1) * n_angles;
    cellmesh.faces(side_face_indices1, 2) = circshift(cellmesh.faces(side_face_indices1, 1), -1);
    cellmesh.faces(side_face_indices1, 3) = (1:n_angles) + i * n_angles;
    side_face_indices2 = side_face_indices1 + n_angles;
    cellmesh.faces(side_face_indices2, 1) = cellmesh.faces(side_face_indices1, 1);
    cellmesh.faces(side_face_indices2, 2) = cellmesh.faces(side_face_indices1, 3);
    cellmesh.faces(side_face_indices2, 3) = circshift(cellmesh.faces(side_face_indices2, 2), 1);
end

% find out which slice to put the nucleus in
% use model for fraction of cell height at which nucleus starts
try
    f_nucbottom = model.nucbottom.stat;
    if f_nucbottom.sigma == 0 
        nucbottom = f_nucbottom.mu;
    else
        nucbottom = ml_rnd(f_nucbottom);
        while tp_pvalue(nucbottom,f_nucbottom) < param.alpha
            nucbottom = ml_rnd(f_nucbottom);
        end
    end
    nucbottomslice = round(cell_height*nucbottom)+1;
    nucbottomslice = min(nucbottomslice, cell_height-nucimgsize(3));
catch
    nucbottomslice = round((cell_height-nucimgsize(3))/2+1);
end


maxsurf = max(nucleus.nucsurf,[],1);
any_points_outside_image = false;

for i = 1:cell_height
    % if current slice should contain some nucleus, calc which slice of
    % nucsurf to use.
    % if not, use the lowest or highest nuclear slice as reference for
    % synthesizing cell shape
    j=max(1,i-nucbottomslice+1);
    j=min(j,nucimgsize(3));
    % use the proper nuclear image slice for reference but clear it
    % after use if above or below where nucleus should be
    
    abovetop = false;
    
    %if we're above the top of the nucleus    
    if i>nucbottomslice+nucimgsize(3)-1
        nucslice = zeros(nucimgsize(1),nucimgsize(2));
        abovetop = true;
    %if we're below the bottom of the nucleus
    elseif i<nucbottomslice 
        nucslice = zeros(nucimgsize(1),nucimgsize(2));
    else
        % make sure that nucleus is inside cell
        % idx=find(nuc2cell_ratio_interp(i,:)>0.95);
        % nuc2cell_ratio_interp(i,idx)=0.95;
        nucslice = nucleus.nucimg(:,:,j);
    end
    
    
    %if we're above the top of the cell, then we enforce that the cell
    %boundary is at least as distant as the nucleus
    if ~abovetop
        nuc2cell_ratio_interp(i,nuc2cell_ratio_interp(i,:)>0.95) = 0.95;
    else
        1;
    end
    
%     [x,y] = pol2cart(Phi,nucleus.nucsurf(j,:)./nuc2cell_ratio_interp(i,:));
    [x,y] = pol2cart(Phi,maxsurf./nuc2cell_ratio_interp(i,:));
    
    x = x + xcenter;
    y = y + ycenter;

    vertices_indices = (1:n_angles) + (i-1) * n_angles;
    cellmesh.vertices(vertices_indices, 1:2) = [x(1:end-1)', y(1:end-1)'];
    cellmesh.vertices(vertices_indices, 3) = i;
    
    [previous_warning_str, previous_warning_id] = lastwarn();
    warning_state = warning;
    warning('off', 'CellOrganizer:ml_interpcurve2img:outOfRange');
    
    %i, min(x), max(x), min(y), max(y)
    sliceimg = ml_interpcurve2img( nucimgsize(1:2), x, y);
    
    [warning_str, warning_id] = lastwarn();
    warning(warning_state);
    if ~strcmp(previous_warning_id, warning_id) && strcmp(warning_id, 'CellOrganizer:ml_interpcurve2img:outOfRange')
        any_points_outside_image = true;
    end

    %D. Sullivan 7/6/13 - doesn't make any sense to get the perim image and
    %then fill it! 
    % if max(nucslice(:))
        % [xn,yn] = pol2cart(Phi,nucleus.nucsurf(j,:));
        % xn = xn + xcenter;
        % yn = yn + ycenter;
        % nucslice = ml_interpcurve2img( nucimgsize(1:2), xn, yn);
       % % figure(1); plot(x,y,xn,yn);
    % end
    % figure(2); imshow(sliceimg); figure(3); imshow(nucslice); pause(0.1);
    cellimg(:,:,i) = imfill(sliceimg,'holes');
    %D. Sullivan 7/6/13 should be no reason to fill again
    % nucimg(:,:,i) = imfill(nucslice,'holes');
    nucimg(:,:,i) = nucslice;
end

if any_points_outside_image
    warning('CellOrganizer:ml_gencellshape3d:outOfRange', 'CellOrganizer: Out of range subscripts when adding pixels to image');
end


result = struct('cellimg', cellimg, 'nucimg', nucimg, 'cellmesh', cellmesh);
