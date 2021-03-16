function [I_transform, param] = tcell_cell_shape_matching(I_seg, I_raw, cellMembrane, model, param)
% The purpose fo the funciton is to map the shape of the template to the shape 
% of the synthesized cell shape. The mapping also uses LDDMM model in CellOrganizer.
% And in mapping we also rescale the image to the cell shape image. 
% 
% After this step, there are the images for cell shape, nuclear shape and protein pattern. 
% And these images are further processed by the functions, the same as other 
% synthesis frameworks in CellOrganizer. 
% 
% Author: Xiongtao Ruan


registration_options = model.t_cell_info.morphing_info.registration_options;

if false 
    registration_options = [];
    registration_options = struct();

    % registration_options.filter_radius = 4;
    % registration_options.filter_radius = 8;
    % registration_options.filter_radius = 16;
    % registration_options.filter_radius = 32;
    % registration_options.filter_radius = ceil(master_script_options.cell_diameter / 4 / 2) * 2;
    registration_options.filter_radius = ceil(master_script_options.cell_diameter / 2 / 2) * 2;
    registration_options.kernel_z_radius = registration_options.filter_radius;
    registration_options.maximum_deformation_per_step = [1, 1, 1] * .5;
    registration_options.window_radius = t_cell_info.template_info.template_options.imx;
    % registration_options.convergence_absolute_error = prod(cropped_size) ./ 2e3;
    % registration_options.convergence_absolute_error = prod(cropped_size) ./ 1e3;
    % registration_options.convergence_absolute_error = prod(cropped_size) ./ 5e2;
    registration_options.convergence_absolute_error = prod(cropped_size) ./ 2e2;
    % registration_options.convergence_absolute_error = prod(cropped_size) ./ 1e2;
    % registration_options.convergence_absolute_error = prod(cropped_size) ./ 1e1;
    registration_options.single_sided = true;
    registration_options.maximum_iterations = max(t_cell_info.template_info.template_dimensions ./ registration_options.maximum_deformation_per_step);

    registration_options.save_intermediates = false;
    registration_options.save_stages = false;

end

% rescale the image if necessary

scale_factor = (sum(cellMembrane(:) > 0) ./ sum(I_seg(:) > 0)) .^ (-1 / 3);

transfrom_matrix = eye(4);
transfrom_matrix(1 : 3, 1 : 3) = eye(3) * scale_factor;


I_cell = affine(cellMembrane, transfrom_matrix([2, 1, 3, 4], [2, 1, 3, 4]), [], false, [], [], 'crop');

I_cell = I_cell > 0.5;

[sy, sx, sz] = size(I_seg);

[I_cell_crop, Bound_box] = crop_image_by_bounding_box(I_cell);
[sy_1, sx_1, sz_1] = size(I_cell_crop);

sy_p = max(sy, sy_1);
sx_p = max(sx, sx_1);
sz_p = max(sz, sz_1);


I_cell_pad = zeros(sy_p, sx_p, sz_p);

I_cell_pad(floor((sy_p - sy_1) / 2) + 1 : floor((sy_p - sy_1) / 2) + sy_1, floor((sx_p - sx_1) / 2) + 1 : floor((sx_p - sx_1) / 2) + sx_1, ...
    floor((sz_p - sz_1) / 2) + 1 : floor((sz_p - sz_1) / 2) + sz_1) = I_cell_crop;

I_seg_pad = zeros(sy_p, sx_p, sz_p);

I_seg_pad(floor((sy_p - sy) / 2) + 1 : floor((sy_p - sy) / 2) + sy, floor((sx_p - sx) / 2) + 1 : floor((sx_p - sx) / 2) + sx, ...
    floor((sz_p - sz) / 2) + 1 : floor((sz_p - sz) / 2) + sz) = I_seg;

I_raw_pad = zeros(sy_p, sx_p, sz_p);

I_raw_pad(floor((sy_p - sy) / 2) + 1 : floor((sy_p - sy) / 2) + sy, floor((sx_p - sx) / 2) + 1 : floor((sx_p - sx) / 2) + sx, ...
    floor((sz_p - sz) / 2) + 1 : floor((sz_p - sz) / 2) + sz) = I_raw;

source = I_seg_pad;

target = I_cell_pad;

r = Greedy3D_lambda_pre(source, target, 1, registration_options)

I_raw_transform = interp3(I_raw_pad, r.source_deformation{1}, r.source_deformation{2}, r.source_deformation{3});

[~, Bbox] = crop_image_by_bounding_box(I_cell_pad);
I_transform_rescale = zeros(size(I_cell));
I_transform_rescale(Bound_box(1, 1) : Bound_box(1, 2), Bound_box(2, 1) : Bound_box(2, 2), Bound_box(3, 1) : Bound_box(3, 2)) =  ...
    I_raw_transform(Bbox(1, 1) : Bbox(1, 2), Bbox(2, 1) : Bbox(2, 2), Bbox(3, 1) : Bbox(3, 2));

transfrom_matrix_invert = eye(4) / transfrom_matrix;

I_transform = affine(I_transform_rescale, transfrom_matrix_invert([2, 1, 3, 4], [2, 1, 3, 4]), [], false, [], [], 'crop');
I_transform = contrast_stretch(I_transform);
end