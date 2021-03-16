function [cell] =  tcell_to_parameter(cell_option, options)
% The script is called in tcell_img2param.m. The script is the main script 
% to parameterize a cell. The input of the function are: cell_options: the 
% option that specific to the cell, such as cell index; options: general
% options that contains the information for the input of the functions. 
% The output is the parameter of the cell. 
%
% Author: Xiongtao Ruan


% extract options for the cell. 
t_cell_info = options.t_cell_info;

% run_key_use_list = t_cell_info.synapse_info.run_key_use_list;
image_name_run_list = t_cell_info.synapse_info.image_name_run_list;
frame_channel_list = t_cell_info.synapse_info.frame_channel_list;
synapse_tracks_run_list = t_cell_info.synapse_info.synapse_tracks_run_list;
relative_time_run_list = t_cell_info.synapse_info.relative_time_run_list;
frame_tracking_list = t_cell_info.synapse_info.frame_tracking_list;
% cell_number_run_list = t_cell_info.synapse_info.cell_number_run_list;  
% run_key_index_list = t_cell_info.synapse_info.run_key_index_list;

i = cell_option.cell_index;

% cell_option.run_key = run_key_use_list{i};
cell_option.image_name = image_name_run_list{i};
cell_option.synapse_annotation = synapse_tracks_run_list(i, :);
cell_option.relative_time = relative_time_run_list(i);
cell_option.frame_channel = frame_channel_list(i);
cell_option.frame_index = frame_tracking_list(i, :);
% cell_option.run_index = run_key_index_list(i);
curr_cell_filename = cell_option.image_name;

% crop the big image to just contain the bounding box of the target cell
[t_cell_info] = tcell_produce_windows(t_cell_info, cell_option);

% segment the cell in the chosen region
[t_cell_info] = tcell_segment_windows(t_cell_info, cell_option);

% rigid alignment of the segmented cell so that the synapse region face upward approximately
[t_cell_info] = tcell_align_segmentations(t_cell_info, cell_option);


% perform alignment adjustment to improve the alignment
if ~t_cell_info.options.use_two_point_synapses && t_cell_info.options.adjust_one_point_alignment
    alignments_adjust_location = [options.temporary_results, '/alignments_adjust/'];
    mkdir(alignments_adjust_location);
    t_cell_info.path_info.alignments_adjust_location = alignments_adjust_location;    
    [t_cell_info] = tcell_adjust_one_point_alignment(t_cell_info, cell_option);
end

% non-rigid alignment of the aligned cell to a half-elipsoid template
[t_cell_info, cell] = tcell_morph_segmentations(t_cell_info, cell_option);    

end