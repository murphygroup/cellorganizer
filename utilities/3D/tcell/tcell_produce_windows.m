function [t_cell_info] = tcell_produce_windows(t_cell_info, options)
% The inputs are t_cell_info, options, the former one contains the general 
% information for the pipeline, and the latter contains the information for 
% the specific cell. The output is  t_cell_info. And the output result relative 
% to the cell is saved in the disk. 
% 
% The idea is to use the coordinates in the annotation file and crop a bounding 
% box that contains the cell with fixed size (71 X 71 X 35). 
% 2016-02-25 xruan: Copied from master_script_produce_windows.m.
% 
% 
% 
% 
% 
% % Dependencies:
% % From the File Exchange: export_fig, pmkmp, toolbox_graph, toolbox_fast_marching, affine, Snake3D, Mesh_voxelisation, dirr, inhull
% % From tebuck: modified version of Jieyue Li's HPA_lib, rasterize_mesh
%
% Author: Xiongtao Ruan



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Produce windows from images corresponding to individual synapses:


% Variables referenced more than once can be copied to local variables:
master_script_options = t_cell_info.options;
regions_location = t_cell_info.path_info.regions_location;
room = t_cell_info.preprocessing_info.room;
image_channels = t_cell_info.preprocessing_info.image_channels;

window_size_2d = t_cell_info.preprocessing_info.window_size_2d;
window_center_2d = t_cell_info.preprocessing_info.window_center_2d;

desired_voxel_size = t_cell_info.synapse_info.desired_voxel_size;
design_voxel_size = t_cell_info.synapse_info.design_voxel_size;

% parameters for the image. 
image_name = options.image_name;
synapse_annotation = options.synapse_annotation;
relative_time = options.relative_time;
frame_channel = options.frame_channel;
frame_index = options.frame_index;
run_index = frame_index(:, 2);

% run_image_size = t_cell_info.synapse_info.run_image_size_list(run_index, :);
run_image_voxel_size = t_cell_info.synapse_info.run_image_voxel_size_list(run_index, :);
current_room = room * desired_voxel_size(1) / run_image_voxel_size(1);
current_z_desired_scale = run_image_voxel_size(3) / desired_voxel_size(3);
% disp(sprintf('Creating window images for run "%s"\n', run_key))

number_windows_written = 0;

current_frame_image_filename = image_name;

% currently only implement one channel case. 
if image_channels == 1
    try
        current_frame_image = robust_read_stack(current_frame_image_filename);
    catch
       fprintf('Unable to read the image: %s\n', current_frame_image_filename);
    end
else
    error('The functionality for multi-channel cases is under development!');
end

% Check if this cell appears in this frame:
current_synapse_annotations = options.synapse_annotation;
if isempty(current_synapse_annotations)
    fprintf('%s annotation is empty', current_synapse_annotations);
    return;
end
% get the synapse center 
current_synapse_center_rounded = round([mean(current_synapse_annotations([1, 3])), mean(current_synapse_annotations([2, 4]))]);

current_relative_time = relative_time;
% generate the image filename based on the run number cell number frame
% number and synapse coordinates
current_window_image_filename = sprintf('run%d_cell%02d_frame%05d_synapse%05d,%05d', frame_index(2), frame_index(3), frame_index(1), current_synapse_center_rounded);
if any(current_synapse_annotations < 0) || any(current_synapse_annotations([1, 3]) > size(current_frame_image, 2)) || any(current_synapse_annotations([2, 4]) > size(current_frame_image, 1))
    fprintf('Synapse annotation coordinate error, skip the preprocess for %s!\n\n', current_window_image_filename);
    save([regions_location, current_window_image_filename], '');
    return;
end
[can_start, final_name, final_exists] = chunk_start_clean(regions_location, current_window_image_filename);
if final_exists
    fprintf('\nCropped window %s already exists\n', current_window_image_filename);
    return;
end
% the process for cropping the image. 
while ~final_exists
    if ~can_start && ~final_exists
        out_status = chunk_lock_clean(regions_location);
        if out_status
            disp('Delete lock files whose jobs are not running!\n')
        end
    end

    fprintf('\nCrop window for frame "%s"\n', current_window_image_filename);        
    current_window_image = imtransform(current_frame_image, maketform('affine', eye(3)), 'bilinear', 'UData', [1, size(current_frame_image, 2)], 'VData', [1, size(current_frame_image, 1)], 'Size', [window_size_2d, size(current_frame_image, 3)], 'XData', current_synapse_center_rounded(1) + current_room * [-1, 1], 'YData', current_synapse_center_rounded(2) + current_room * [-1, 1], 'FillValues', nan);
    % Rescale 
    current_window_image = image_resize_nd(current_window_image, [1, 1, current_z_desired_scale], 'bilinear', true);
    % Transform annotations for alignment convenience:
    % current_synapse_center_rounded_transformed = window_center_2d;
    % N x 2:
    current_synapse_annotations_transformed = reshape(current_synapse_annotations, 2, []).';
    % Translate and scale:
    current_synapse_annotations_transform = struct('current_synapse_center_rounded', current_synapse_center_rounded, 'current_room', current_room, 'room', room);
    current_synapse_annotations_transform_function = @(given_annotations)(given_annotations - repmat(current_synapse_annotations_transform.current_synapse_center_rounded, [size(given_annotations, 1), 1]) + current_synapse_annotations_transform.current_room) .* (current_synapse_annotations_transform.room / current_synapse_annotations_transform.current_room) + 1;
    current_synapse_annotations_transformed = current_synapse_annotations_transform_function(current_synapse_annotations_transformed);
    current_synapse_center_rounded_transformed = current_synapse_annotations_transformed(1, :);

    current_window_image = extrapolate_nans(current_window_image);

    save([regions_location, current_window_image_filename, '.mat'], 'current_window_image', 'current_synapse_center_rounded', 'current_synapse_annotations', 'current_synapse_annotations_transform', 'current_synapse_center_rounded_transformed', 'current_synapse_annotations_transformed', 'current_relative_time')

    number_windows_written = number_windows_written + 1;
    chunk_finish(regions_location, current_window_image_filename);  
    final_exists = true;
end
% fsprintf('  Created %d window images\n', number_windows_written)
t_cell_info.produce_windows_info = struct();
t_cell_info.produce_windows_info.room = room; 


end  

