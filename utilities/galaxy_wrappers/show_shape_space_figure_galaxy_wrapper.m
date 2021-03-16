function answer = show_shape_space_figure_galaxy_wrapper( model_filename, options )
%SHOW_SHAPE_SPACE_FIGURE_GALAXY_WRAPPER Helper function that reads a model and saves
%an RGB png file of the shape space.

% Ivan E. Cao-Berg (icaoberg@cmu.edu)
%
% Copyright (C) 2016 Murphy Lab
% Carnegie Mellon University
%
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

% specify model name and load model
if ~exist( model_filename )
    answer = false;
    warning(['Model filename ' model_filename ' was not found on disk.']);
    return
else
    load(model_filename);
end

%%%%%
% demo2d05:pca
% demo2d04:diffeomorphic
% demo3D51:spharm
%%%%%

if strcmp(model.nuclearShapeModel.type,'diffeomorphic')
    nimgs = size(model.cellShapeModel.positions,1);
    if options.nlabels > nimgs
        answer = false;
        warning(['The number of labels you requested (' ...
            num2str(options.nlabels) ') is greater than the number of images ' ...
            '(' num2str(nimgs) ') available. Setting the number of ' ...
            'labels to the number of available images.']);
        options.nlabels = nimgs;
    end
    labels = reshape(repmat([1:options.nlabels],ceil(nimgs/options.nlabels),1),[],1);
    options.plot_dims = [1,2];
    % show shape space by calling the function
    f = figure('visible','off');
    showShapeSpaceFigure( model, labels, options ); %% this is diffeomorphic
elseif strcmp(model.nuclearShapeModel.type,'pca')
    % show shape space by calling the function
    f = figure('visible','off');
    showPCAShapeSpaceFigure(model);
elseif strcmp(model.nuclearShapeModel.type,'spharm_rpdm')
    % show shape space by calling the function
    f = figure('visible','off');
    show_SPHARM_RPDM_Shape_Space_Figure(model);
else
    warning('This is not a PCA, diffeomorphic, or shparm model, cannot show shape space figure.');
    return
end

saveas( f, 'show_shape_space.png', 'png' );

end%show_shape_space_figure_galaxy_wrapper