function [t_cell_info] = tcell_get_morphing_info(t_cell_info)
    master_script_options = t_cell_info.options;
    regions_location = t_cell_info.path_info.regions_location;
    segmentations_filtered_location = t_cell_info.path_info.segmentations_filtered_location;
    alignments_location = t_cell_info.path_info.alignments_location;
    deformations_location = t_cell_info.path_info.deformations_location;

    template_centroid = t_cell_info.template_info.template_centroid;
    template_synapse = t_cell_info.template_info.template_synapse;
    template_image = t_cell_info.template_info.template_image;

    window_size_2d = t_cell_info.preprocessing_info.window_size_2d;
    window_center_2d = t_cell_info.preprocessing_info.window_center_2d;
    cropped_size = t_cell_info.preprocessing_info.cropped_size;

    % Set options for T cell image registration:
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
    % registration_options.save_intermediate_images = true;
    % registration_options.always_save_full_intermediates = true;
    % registration_options;
    % beep, keyboard
    cell_param = struct();

    t_cell_info.morphing_info = struct();
    t_cell_info.morphing_info.registration_options = registration_options;
  
end