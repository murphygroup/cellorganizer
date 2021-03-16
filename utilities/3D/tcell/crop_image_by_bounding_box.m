function [I_crop, Bound_box] = crop_image_by_bounding_box(I)

I_y = squeeze(sum(sum(I, 2), 3));

I_x = squeeze(sum(sum(I, 1), 3));

I_z = squeeze(sum(sum(I, 1), 2));

y1 = find(I_y > 0, 1, 'first');
y2 = find(I_y > 0, 1, 'last');

x1 = find(I_x > 0, 1, 'first');
x2 = find(I_x > 0, 1, 'last');

z1 = find(I_z > 0, 1, 'first');
z2 = find(I_z > 0, 1, 'last');

I_crop = I(y1 : y2, x1 : x2, z1 : z2);

Bound_box = [y1, y2; x1, x2; z1, z2];

end