function [t_cell_info, cell_param] = tcell_morph_segmentations(t_cell_info, options)
% The inputs are t_cell_info, options, the former one contains the general 
% information for the pipeline, and the latter contains the information for
% the specific cell. The output is  t_cell_info. And the output result relative 
% to the cell is saved in the disk.
% 
% The idea is to use the aligned image in the disk, and morph the cell to the
% template. After rigid alignment, the cell is approximately similar to the 
% template. To further exploring the relationship of the cell and the template, 
% we use a registration method LDDMM model which can map the one image to another 
% image by wrapping the image. After that, we can project the voxels in the 
% aligned image to the voxels in the template, so that different cells can be comparable. 
% 
% After this step, we parameterize the image with the standardized shape. And 
% the next step is to build the model using the parameters of all cells.%
% 
% 2016-02-29 xruan: Copied from master_script_morph_segmentations.m.
% 2017-05-28 xruan: add verbose functionality
% 2018-08-21 xruan: add option for alignment adjustment
% 
% 
% 
% % Dependencies:
% % From the File Exchange: export_fig, pmkmp, toolbox_graph, toolbox_fast_marching, affine, Snake3D, Mesh_voxelisation, dirr, inhull
% % From tebuck: modified version of Jieyue Li's HPA_lib, rasterize_mesh
%
% Author: Xiongtao Ruan
  
  
  
% Variables referenced more than once can be copied to local variables:
master_script_options = t_cell_info.options;
regions_location = t_cell_info.path_info.regions_location;
segmentations_filtered_location = t_cell_info.path_info.segmentations_filtered_location;
alignments_location = t_cell_info.path_info.alignments_location;
deformations_location = t_cell_info.path_info.deformations_location;
illustration_location = t_cell_info.path_info.illustration_location;

template_centroid = t_cell_info.template_info.template_centroid;
template_synapse = t_cell_info.template_info.template_synapse;
template_image = t_cell_info.template_info.template_image;

window_size_2d = t_cell_info.preprocessing_info.window_size_2d;
window_center_2d = t_cell_info.preprocessing_info.window_center_2d;
cropped_size = t_cell_info.preprocessing_info.cropped_size;

segmentation_rasterization_function = t_cell_info.segmentation_info.segmentation_rasterization_function;

template_volume = t_cell_info.alignment_info.template_volume;
segmentation_volume_threshold = t_cell_info.alignment_info.segmentation_volume_threshold;

registration_options = t_cell_info.morphing_info.registration_options;

% debug_print_deformation_skip_reasons = false;
debug_print_deformation_skip_reasons = true;

image_name = options.image_name;
synapse_annotation = options.synapse_annotation;
relative_time = options.relative_time;
frame_channel = options.frame_channel;
frame_index = options.frame_index;
run_index = frame_index(:, 2);  


% fprintf_diary_cycle_function('Segmenting window images for run "%s"\n', run_key)
printed_run_work_message = false;

% Go through frames of each cell:
% for current_cell_index = reshape(find(any(cell_index_frame_indices_uncomputed{run_index}, 2)), 1, [])
current_synapse_annotations = options.synapse_annotation;
if isempty(current_synapse_annotations)
    fprintf('%s annotation is empty\n', current_synapse_annotations);
    return;
end
current_synapse_center_rounded = round([mean(current_synapse_annotations([1, 3])), mean(current_synapse_annotations([2, 4]))]);
current_relative_time = relative_time;
cell_deformation_filename = [sprintf('run%d_cell%02d_frame%05d_synapse%05d,%05d', frame_index(2), frame_index(3), frame_index(1), current_synapse_center_rounded)];
cell_alignment_filename = [alignments_location, cell_deformation_filename];
cell_segmentation_filename = [segmentations_filtered_location, cell_deformation_filename];
cell_region_filename = [regions_location, cell_deformation_filename];

cell_param.deformation_filename = cell_deformation_filename;  

[can_start, final_name, final_exists] = chunk_start_clean(deformations_location, cell_deformation_filename);
if final_exists
    fprintf('Morphed cell %s already exists\n', cell_segmentation_filename);
    return;
end

while ~final_exists
    if ~can_start && ~final_exists
        out_status = chunk_lock_clean(deformations_location);
        if out_status
            disp('Delete lock files whose jobs are not running!\n')
        end
    end
    try 
        a = load([cell_segmentation_filename, '.mat']);
        % Feb. 14, 2016 xruan disable it
        % segmentation_structure = a.segmentation_structure;
        window_number_slices = a.window_number_slices;

        % 08/21/2018 add option for alignment adjustment
        if t_cell_info.options.adjust_one_point_alignment
            cell_alignment_adjust_filename = [t_cell_info.path_info.alignments_adjust_location, cell_deformation_filename];
            a = load([cell_alignment_adjust_filename, '.mat']);
        else
            a = load([cell_alignment_filename, '.mat']);                
        end
        current_transform = a.current_transform;
        cropped_segmentation_image = double(a.cropped_segmentation_image >= 0.5);  
        % cropped_segmentation_image = a.cropped_segmentation_image;
        cropped_raw_image = a.cropped_raw_image;            
    catch
        temp = struct; save(final_name, '-struct', 'temp'); chunk_finish(deformations_location, cell_deformation_filename); return;
    end
    % catch raised_error
    % getReport(raised_error, 'extended')
    % end
    if isempty(fieldnames(a))
        temp = struct; save(final_name, '-struct', 'temp'); chunk_finish(deformations_location, cell_deformation_filename);
        if debug_print_deformation_skip_reasons, fprintf('Skipping morph %s cell %d frame %d because alignment empty\n', run_key, cell_number, frame_index), end
        return;
    end
    % xruan 08/14/2015

    current_volume = sum(cropped_segmentation_image(:));
    is_segmentation_low_volume = current_volume <= segmentation_volume_threshold;

    if is_segmentation_low_volume
    % if debug_print_deformation_skip_reasons, fprintf_diary_cycle_function('Skipping morph %s cell %d because cropped_segmentation_image blank\n', run_key, current_cell_index), end
    if debug_print_deformation_skip_reasons, fprintf('Skipping morph %s cell %d frame %d because cropped_segmentation_image blank\n', run_key, cell_number, frame_index), end
        temp = struct; save(final_name, '-struct', 'temp'); chunk_finish(deformations_location, cell_deformation_filename);
        return
    end

    % xruan 02/24/2016 
    fprintf('Morphing aligned segmentations for cell "%s"\n', cell_deformation_filename)

    a = load([cell_region_filename, '.mat']);
    current_window_image = a.current_window_image;

    source = cropped_segmentation_image;
    target = template_image;

    % using diffeomorphic model for the morphing
    r = Greedy3D_lambda_pre(...
    source, target, 1, registration_options)

    cropped_segmentation_image_deformed = interp3(cropped_segmentation_image, r.source_deformation{1}, r.source_deformation{2}, r.source_deformation{3});
    cropped_raw_image_deformed = interp3(cropped_raw_image, r.source_deformation{1}, r.source_deformation{2}, r.source_deformation{3});
    save(final_name, 'cropped_segmentation_image_deformed', 'cropped_raw_image_deformed', 'registration_options', 'r', 'current_window_image')
    chunk_finish(deformations_location, cell_deformation_filename);

    final_exists = true;     
end

% add verbose for the pipeline
if master_script_options.verbose
    close all
    a = load(cell_alignment_filename);
    cropped_raw_image = a.cropped_raw_image;
    cropped_segmentation_image = a.cropped_segmentation_image;

    image_to_show = [];
    image_to_show = [image_to_show; contrast_stretch(cropped_raw_image)];
    image_to_show = [image_to_show; cropped_segmentation_image];
    image_to_show = [image_to_show; contrast_stretch(cropped_raw_image .* (cropped_segmentation_image >= .5) + min(cropped_raw_image(cropped_segmentation_image >= .5)) .* (cropped_segmentation_image < .5))];
    image_to_show = [image_to_show; contrast_stretch(cropped_raw_image_deformed .* (cropped_segmentation_image_deformed >= .5) + min(cropped_raw_image_deformed(cropped_segmentation_image_deformed >= .5)) .* (cropped_segmentation_image_deformed < .5))];
    image_to_show = reshape_2d(image_to_show, 1);
    cell_illustration_filename =[illustration_location, cell_deformation_filename];
    image_filename = cell_illustration_filename;
    imwrite(image_to_show, [image_filename, '.png'])
    fprintf('Wrote illustration "%s"\n\n', cell_deformation_filename)
    figure, imshow(image_to_show, []), pause(1)
end
    
    
end  

