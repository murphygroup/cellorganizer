function [I_geometry1, I_geometry2] = ellipsoid_geometry( model )
% Use a plane to get a half-ellipsoid geometry on model
%
% Input Variables from Model
% --------------------------
% * image_size: size of the image
% * ellipsoid_param: a, b, c for the 3D ellipsoid
% * centroid: centroid of the ellipsoid, default uses the center of the image
% * line_param: linear function that cuts the ellipsoid
%
% Output Variables
% ----------------
% * I_geometry1: output image for the generated geometry, represented as
% binary image
% * I_geometry2: same as I_geometry1
%
% Author: Xiongtao Ruan

image_size =        model.cellShapeModel.csgo.image_size;
ellipsoid_param =   model.cellShapeModel.csgo.shape_param;
centroid =          model.cellShapeModel.csgo.centroid;
plane_param =       model.cellShapeModel.csgo.plane_param;

if nargin < 1
    image_size = [71, 71, 35];
    ellipsoid_param = [26, 26, 14];
    centroid = [36, 36, 18];
    plane_param = [0, -1, 0, -36];
end

% I = zeros(image_size);

[I_y, I_x, I_z] = ndgrid(1:image_size(1), 1:image_size(2), 1:image_size(3));

% parameter for the ellipsoid
% (x - x1)^2 / a^2 + (y - y1)^2 / b^2 + (z - z1)^2 / c^2 <= 1
a = ellipsoid_param(1);
b = ellipsoid_param(2);
c = ellipsoid_param(3);

x1 = centroid(2);
y1 = centroid(1);
z1 = centroid(3);

% parameter for the plane
% k1 * x + k2 * y + k3 * z <= d
k1 = plane_param(1);
k2 = plane_param(2);
k3 = plane_param(3);
d = plane_param(4);


I_ellipsoid = (I_x - x1) .^ 2 / a^2 + (I_y - y1) .^ 2 / b^2 + (I_z - z1) .^ 2 / c^2;

I_ellipsoid_bw = I_ellipsoid <= 1;

I_plane = k1 * I_x + k2 * I_y + k3 * I_z;

I_plane_bw = I_plane <= d;

I_geometry1 = I_ellipsoid_bw .* I_plane_bw;
I_geometry2 = I_geometry1;

if ~false
    figure, imshow(reshape(I_geometry1, size(I_geometry1, 1), []), [])
end

end