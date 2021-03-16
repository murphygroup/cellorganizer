function [I_region] = Tcell_user_defined_synapse_region(template_image, radius, thickness, top_start_ind)
% The function defines a cylinder with specific raidus and height
% (thickness), at the top of the half-ellipse template with the same axis.
% top_start_ind controls the start pixel from the template, that is, the
% number of slices at xz plance not included. 

% Author: Xiongtao Ruan
% Date: 11/04/2018


if nargin < 1
    a = load('standardized_voxels.mat');
    template_image = a.model.template_info.template_image;
    radius = 9;
    thickness = 4;
    top_start_ind = 2;
end

x_range = find(sum(sum(template_image, 1), 3) > 0);
y_range = find(sum(sum(template_image, 2), 3) > 0);
z_range = find(sum(sum(template_image, 1), 2) > 0);

synapse_center = [mean(x_range), y_range(1), mean(z_range)];

region_center = synapse_center + [0, top_start_ind - 1 + thickness / 2, 0];

[X, Y, Z] = meshgrid(1 : size(template_image, 2), 1 : size(template_image, 1), 1 : size(template_image, 3));

% define the region by the function of a cylinder at center (a, b, c) with
% radius r and height h, with axis along y-axis, (x - a)^2 + (z - c)^ 2 <=
% r^2 and |y - b| <= h / 2
I_1 = (X - region_center(1)) .^ 2 + (Z - region_center(3)) .^ 2 <= radius ^ 2;
I_2 = abs(Y - region_center(2)) <= thickness / 2;

I_region = I_1 .* I_2;

I_region = I_region .* template_image > 0;

end