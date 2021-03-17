function answer = demo3D27()
% demo3D27
%
% This demo performs a regression between two sets of related shapes (i.e.
% predicts cell  shape from nuclear shape) and displays the residuals as in
% Figure 2 of Johnson et al 2015.
%
% Input 
% -----
% * models hela_cell_10_15_15.mat and hela_nuc_10_15_15.mat
%
% Output
% ------
% * shape space figure

% Copyright (C) 2016-2017 Murphy Lab
% Computational Biology Department
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
% DO NOT MODIFY THIS BLOCK
if ~isdeployed()
  current_path = which(mfilename);
  [current_path, filename, extension] = fileparts( current_path );
  cd(current_path);
end

disp( 'demo3D27' );
warning('This demo is deprecated. The demo will be removed in future versions of CellOrganizer');
disp( 'The estimated running time is 2 hours. Please wait.' );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT MODIFY THIS BLOCK
current_path = which(mfilename);
[current_path, filename, extension] = fileparts( current_path );
cd(current_path);

modelpath = which(mfilename);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modelCell = load('../../../models/3D/diffeomorphic/hela_cell_10_15_15.mat');
modelNuc = load('../../../models/3D/diffeomorphic/hela_nuc_10_15_15.mat');

distances_incomplete_cell = modelCell.model.cellShapeModel.distances_incomplete;
distances_incomplete_nuc = modelNuc.model.cellShapeModel.distances_incomplete;

[ keepinds_cell ] = get_complete_distance_matrix( distances_incomplete_cell);
[ keepinds_nuc ] = get_complete_distance_matrix( distances_incomplete_nuc);

keepinds = intersect(keepinds_cell, keepinds_nuc);

cellpos = modelCell.model.cellShapeModel.positions(keepinds, :);
nucpos = modelNuc.model.cellShapeModel.positions(keepinds, :);

ndat = size(cellpos,1);

c2ndir = [pwd filesep 'c2n'];
n2cdir = [pwd filesep 'n2c'];

c2n = kernel_pred_iter(cellpos, nucpos, 1:ndat, 0, c2ndir);
n2c = kernel_pred_iter(nucpos, cellpos, 1:ndat, 0, n2cdir);

[cell2nuc_dist_pval, cell2nuc_pred_norm_err] = pred_pval(nucpos, c2n.y_pred);
[nuc2cell_dist_pval, nuc2cell_pred_norm_err] = pred_pval(cellpos, n2c.y_pred);

imfunc_cell = modelCell.model.cellShapeModel.imfunc;
imfunc_nuc = modelNuc.model.cellShapeModel.imfunc;

param.subsize = 150;

pmax = max([cell2nuc_dist_pval; nuc2cell_dist_pval]);

boundpad = 0.25;

positions_nuc_norm = zeros(size(nucpos));
positions_nuc_norm(:) = zscore(nucpos(:));

positions_cell_norm = zeros(size(cellpos));
positions_cell_norm(:) = zscore(cellpos(:));

colors = nan(ndat,3);
colors(keepinds,:) = grayrange2rgb(cell2nuc_dist_pval, 0,pmax);

model = modelCell.model;
model.cellShapeModel.positions =  nan(size(model.cellShapeModel.positions,1), size(positions_nuc_norm,2));
model.cellShapeModel.positions(keepinds,:) = positions_nuc_norm;

param.imfunc = @(x) (im2projection_RGB(permute(cropImg(imfunc_cell(x),0),[2,1,3]),struct('scale_inten', 2, 'justz', true, 'cm', @(y) colormap(colors(x,:)))));

mkdir img;

h(1) = figure;
showShapeSpaceFigure(model,[], param);

minpos = min(model.cellShapeModel.positions,[],1);
maxpos = max(model.cellShapeModel.positions,[],1);
axis([minpos(1)-boundpad maxpos(1)+boundpad minpos(2)-boundpad maxpos(2) + boundpad])

param.imfunc = @(x) (im2projection_RGB(permute(cropImg(imfunc_nuc(x),0),[2,1,3]),struct('scale_inten', 2, 'justz', true, 'cm', @(y) colormap(colors(x,:)))));
saveas( gcf, './img/output01.png' );

h(2) = figure;
showShapeSpaceFigure(model,[], param);

axis([minpos(1)-boundpad maxpos(1)+boundpad minpos(2)-boundpad maxpos(2) + boundpad]);

colors = nan(ndat,3);
colors(keepinds,:) = grayrange2rgb(nuc2cell_dist_pval, 0,pmax);

model.cellShapeModel.positions = nan(size(model.cellShapeModel.positions,1), size(positions_cell_norm,2));
model.cellShapeModel.positions(keepinds,:) = positions_cell_norm;

param.imfunc = @(x) (im2projection_RGB(permute(cropImg(imfunc_cell(x),0),[2,1,3]),struct('scale_inten', 2, 'justz', true, 'cm', @(y) colormap(colors(x,:)))));
saveas( gcf, './img/output02.png' );

h(3) = figure;
showShapeSpaceFigure(model,[], param);

minpos = min(model.cellShapeModel.positions,[],1);
maxpos = max(model.cellShapeModel.positions,[],1);
axis([minpos(1)-boundpad maxpos(1)+boundpad minpos(2)-boundpad maxpos(2) + boundpad])

param.imfunc = @(x) (im2projection_RGB(permute(cropImg(imfunc_nuc(x),0),[2,1,3]),struct('scale_inten', 2, 'justz', true, 'cm', @(y) colormap(colors(x,:)))));
saveas( gcf, './img/output03.png' );

h(4) = figure;
showShapeSpaceFigure(model,[], param);
axis([minpos(1)-boundpad maxpos(1)+boundpad minpos(2)-boundpad maxpos(2) + boundpad]);
saveas( gcf, './img/output04.png' );
answer = true;
end%demo3D27