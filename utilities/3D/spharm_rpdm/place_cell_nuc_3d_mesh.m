function [ h ] = place_cell_nuc_3d_mesh( fvec, location, param )
%places an image centered at the middle most pixel at the coordinate
%specified by the 1d vector location
%
%by default, all zero pixels are set to tranparent
%
%returns the handle for the placed image

% xruan 04/29/2018 

if ~exist('param', 'var')
    param = [];
end

param = ml_initparam(param, struct( ...
        'subsize', 400, ...
        'alpha', [0,0,0], ...
        'colormap_func', @jet, ...
        'viewangle', [45,45] ...
        ));

cell_vertices = [];
nuc_vertices = [];
if isfield(fvec, 'cell_vertices') 
    cell_vertices = fvec.cell_vertices;
    cell_faces = fvec.cell_faces;
end

if isfield(fvec, 'nuc_vertices')
    nuc_vertices = fvec.nuc_vertices;
    nuc_faces = fvec.nuc_faces;
end

% get the center of mesh
if isempty(cell_vertices)
    mesh_center = mean(nuc_vertices);
else
    mesh_center = mean(cell_vertices);
end

if ~isempty(cell_vertices)
    cell_vertices = cell_vertices - mesh_center;
    cell_vertices = cell_vertices + location;
end

if ~isempty(nuc_vertices)
    nuc_vertices = nuc_vertices - mesh_center;
    nuc_vertices = nuc_vertices + location;
end

if ~isempty(cell_vertices)
    if strcmpi(class(param.colormap), 'function_handle')
        h = patch('vertices', cell_vertices, 'faces', cell_faces, 'FaceVertexCData',param.colormap(size(cell_vertices,1)),'FaceColor','interp', 'EdgeColor', 'None');
    else
        h = patch('vertices', cell_vertices, 'faces', cell_faces, 'FaceColor', param.colormap, 'EdgeColor', 'None');
    end
    if ~isempty(nuc_vertices)
        alpha(0.5);
        hold on
        h = patch('vertices', nuc_vertices, 'faces', nuc_faces, 'FaceColor', '#AAAAAA', 'EdgeColor', 'None');
    end
elseif ~isempty(nuc_vertices)
    h = patch('vertices', nuc_vertices, 'faces', nuc_faces, 'FaceVertexCData',param.colormap(size(nuc_vertices,1)),'FaceColor','interp', 'EdgeColor', 'None');
end

view(param.viewangle);
    
% alphamap = zeros(size(image));
% for i = 1:length(param.alpha)
%     alphamap(:,:,i) = param.alpha(i);
% end

% set(h, 'AlphaDataMapping', 'none')
% set(h, 'AlphaData', ~all(image == alphamap,3))


end

