function [ h ] = place_3d_mesh( fvec, location, param )
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
        'colormap_func', @jet ...
        ));

%scale such maximum is 1    
vertices = fvec.vertices;
faces = fvec.faces;
vertices = vertices - mean(vertices);
vertices = vertices + location;

% h = patch('vertices', vertices, 'faces', faces, 'FaceVertexCData',jet(size(vertices,1)),'FaceColor','interp', 'EdgeColor', 'None');
h = patch('vertices', vertices, 'faces', faces, 'FaceVertexCData',param.colormap_func(size(vertices,1)),'FaceColor','interp', 'EdgeColor', 'None');
view([45, 45]);

% alphamap = zeros(size(image));
% for i = 1:length(param.alpha)
%     alphamap(:,:,i) = param.alpha(i);
% end

% set(h, 'AlphaDataMapping', 'none')
% set(h, 'AlphaData', ~all(image == alphamap,3))


end

