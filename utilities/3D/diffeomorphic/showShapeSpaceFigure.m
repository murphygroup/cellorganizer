function [param] = showShapeSpaceFigure(model, labels, param)
% SHOWSHAPESPACEFIGURE 

% Greg Johnson
%
% Copyright (C) 2015-2016 Murphy Lab
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

%need to add:
%   traces - grj 11/15/14
%   scaling option

% $Revision: 0.1
% 20151226 icaoberg Updated error message and included number of objects as
% title

% $Revision: 0.2
% 20160625 rfmurphy added traces (lines between cells of the same label)

if ~exist('param', 'var')
    param = [];
end

nimgs = size(model.cellShapeModel.positions,1);

if size(model.cellShapeModel.positions,2) >= 3
    [~, synthorder] = sort(model.cellShapeModel.positions(:,3), 'Ascend');
else
    synthorder = 1:nimgs;
end

param = ml_initparam(param, struct( ...
    'skipmissing', true, ...
    'rebuild_pos', false, ...
    'plot_dims', [1,2], ...
    'cm', @jet, ...
    'traces', [], ...
    'embedfctn', @cmdscale, ...
    'subsize', 400, ...
    'synthorder', synthorder, ...
    'rotate', 0 ...
    ));
cellnums = 1:nimgs;

%logic: Check for model.dataset.labels 
%       => if not there check for non empty labels in options structure
%       => if not there set it 1 for an array of the length of images
if isfield(model, 'dataset') && isfield(model.dataset, 'labels') && ~isempty(model.dataset.labels)
    [ulabels, ~, labels] = unique(model.dataset.labels);
elseif exist('labels', 'var') && ~isempty(labels)
    [ulabels, ~, labels] = unique(labels);
else
    labels  = ones(1, nimgs); %unique(1:nimgs);
    ulabels = [1]; 
end

skipmissing = param.skipmissing;
cm = param.cm;
plot_dims = param.plot_dims;
rebuild_pos = param.rebuild_pos;

ndims = 2;

if rebuild_pos
    if skipmissing
        %remove rows fully populated with nans
        d = model.cellShapeModel.distances_incomplete;
        d_nan = d;
        d_nan(logical(eye(size(d)))) = nan;
        
        rminds = all(isnan(d_nan),1);
        d(rminds,:) = [];
        d(:,rminds) = [];
        
        keepinds = find(~rminds);
        
        labels = labels(keepinds);
        cellnums = cellnums(keepinds);
        
        %remove any other entries with nans
        rminds = any(isnan(triu(d)),1);
        d(rminds,:) = [];
        d(:,rminds) = [];
        
        keepinds = ~rminds;
        
        labels = labels(keepinds);
        cellnums = cellnums(keepinds);
    else
        d = model.cellShapeModel.distances;
    end
    [embed_pos, eig] = param.embedfctn(d);
else
    eig = [];
    keepinds = ~any(isnan(model.cellShapeModel.positions),2);
    %     labels = labels(keepinds);
    cellnums = cellnums(keepinds);
    %     embed_pos = model.cellShapeModel.positions(keepinds,:);
end


if strcmpi(class(cm), 'function_handle')
    colors = cm(length(ulabels))*0.8;
else
    colors = cm;
end

param.colors = colors;

% if exist('traces', 'var')
%     traces = traces(all(ismember(traces, keepinds),2),:);
%
%     tracemap = ones(max(keepinds), 1) * -1;
%     tracemap(keepinds) = 1:length(keepinds);
%
%
%     shapeInserter = vision.ShapeInserter('Shape','Lines','BorderColor','Custom', 'CustomBorderColor', uint8([125 125 125]));
%
%    lines = int32([locs(tracemap(traces(:,1)),2:-1:1) locs(tracemap(traces(:,2)),2:-1:1)]);
%
%    lines = lines + int32(repmat(imsize./4, [size(traces,1) 2]));
% %         line = int32([x(i),y(i),x(j), y(j)]);
%     img = step(shapeInserter, img, lines);
% end

set(gcf, 'color', 'w')
hold on


param = ml_initparam(param, struct( ...
    'imfunc', @(x) imrotate(img2flat(model.cellShapeModel.imfunc(x), colors(labels(x),:)), param.rotate) ...
    ));

if (param.traces)
    for cellind = 1:length(cellnums)-1
        if labels(cellind)==labels(cellind+1)
            xpos(1) = model.cellShapeModel.positions(cellind,param.plot_dims(1));
            ypos(1) = model.cellShapeModel.positions(cellind,param.plot_dims(2));
            xpos(2) = model.cellShapeModel.positions(cellind+1,param.plot_dims(1));
            ypos(2) = model.cellShapeModel.positions(cellind+1,param.plot_dims(2));
            plot(xpos,ypos,'k-');
        end
    end
end

param.synthorder = param.synthorder(ismember(param.synthorder, cellnums));
for imnum = 1:length(param.synthorder)
    
    cellind = param.synthorder(imnum);
    try
        cellimg = param.imfunc(cellind);
        
        if isempty(cellimg) %|| any(cellimg(:) .* cellbounds(:))
            continue;
        end
        
        %center the image on the point where it belongs
        place_image(cellimg, model.cellShapeModel.positions(cellind,param.plot_dims), param);
        hold on
    catch err
        disp(['Unable to plot image: ' num2str(imnum)] );
        getReport( err )
    end
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