function [ h ] = img2vol( img , scale, alpha)
%IMG2VOL Summary of this function goes here
%   Detailed explanation goes here

if ndims(img) > 3
    warning('Unsupported dimension');
    return;
end

if ~exist('scale', 'var') || isempty(scale)
    scale = [1, 1, 1];
end

img = double(img);
if ~exist('alpha', 'var') || isempty(alpha)
    alpha = img./max(img(:)) / 5;
end
%     img = img(:,:,end:-1:1);

    

    model.cdata = img;
    model.xdata = [0, size(img,2)* scale(2)];
    model.ydata = [0, size(img,1)* scale(1)];
    model.zdata = [0, size(img,3)* scale(3)];
    
    
    model.texture = '3D';
   
    model.parent = [];
    model.tag = '';    
    alphamap = ones(size(img)) .* alpha;
    alphamap(img == 0) = 0;
    
    model.alpha = alphamap;
    
    h = vol3d(model);
    axis tight;
    axis equal;
    set(gcf, 'Color', 'w')
    axis tight;
    h = h.parent;
end

