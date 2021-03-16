function [cell_images, options] = tcell_setup_options(options)
% xruan 2016-03-06 xruan: copied and modified from master_script.m in t
% cell pipeline
% it is called by setup_data.m in img2model.m. The function of the script 
% is to set up the basic options as input for the further steps, such as, 
% making directories in the temp directories, getting template information,
% getting parameters for segmentation, get parameters for morphing and so on. 
% 
% Author: Xiongtao Ruan
% 
% 2017-05-28 xruan: add illustration location 
% 2018-08-20 xruan: add function for automatic synapse detection. 
  
  
  
warning('off', 'MATLAB:MKDIR:DirectoryExists');
warning('off', 'Images:initSize:adjustingMag');
warning('off', 'MATLAB:xlsread:Mode');
warning('off', 'MATLAB:tifftagsread:rational:badTagValueDivisionByZero');

default_options = tcell_get_default_options;
  
if ~exist('options', 'var')
    options = default_options;
else
    options_prior_to_processing = options;
    options = tcell_process_options(options);
end
options
named_options_set = options.named_options_set;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start storing important information in one big structure that can be passed around (instead of using the base namespace or stuff like that):

t_cell_info = struct();
t_cell_info.options = options;

function [result_strings] = condition_name_abbreviation_function(given_strings)
    result_strings = given_strings;
    for abbreviation_index = 1:size(options.condition_name_abbreviations, 1)
        result_strings = strrep(result_strings, options.condition_name_abbreviations{abbreviation_index, 1}, options.condition_name_abbreviations{abbreviation_index, 2});
    end
end

t_cell_info.abbreviation_info.condition_name_abbreviation_function = @condition_name_abbreviation_function;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Paths for saving intermediate and final results:

% Build up a prefix that allows maximal reuse of common results but separation of chains of results that use a different method in an earlier link:
accumulated_directory_options_prefix = '';

should_use_short_directory_names = false;
%should_use_short_directory_names = true

temp_location = [options.tmp_results_location, accumulated_directory_options_prefix, 'temp/'];
if should_use_short_directory_names
    cell_positions_location = [options.tmp_results_location, accumulated_directory_options_prefix, 'cellpos/'];
    regions_location = [options.tmp_results_location, accumulated_directory_options_prefix, 'reg/'];
    segmentations_location = [options.tmp_results_location, accumulated_directory_options_prefix, 'seg/'];
else
    cell_positions_location = [options.tmp_results_location, accumulated_directory_options_prefix, 'cell-positions/'];
    regions_location = [options.tmp_results_location, accumulated_directory_options_prefix, 'regions/'];
    segmentations_location = [options.tmp_results_location, accumulated_directory_options_prefix, 'segmentations/'];
end
if options.use_segmentation_filtering
    segmentation_filtering_interactive_prefix = '';
    if ~options.segmentation_filtering_interactive
        if should_use_short_directory_names
            segmentation_filtering_interactive_prefix = '-man';
        else
            segmentation_filtering_interactive_prefix = '-manual';
        end
    end
    if should_use_short_directory_names
        segmentation_evaluations_location = [options.tmp_results_location, accumulated_directory_options_prefix, 'segeval', segmentation_filtering_interactive_prefix, '/'];
        segmentations_filtered_location = [options.tmp_results_location, accumulated_directory_options_prefix, 'segfiltd', segmentation_filtering_interactive_prefix, '/'];
        accumulated_directory_options_prefix = [accumulated_directory_options_prefix, 'segfiltg', segmentation_filtering_interactive_prefix, '_'];
    else
        segmentation_evaluations_location = [options.tmp_results_location, accumulated_directory_options_prefix, 'segmentation-evaluations', segmentation_filtering_interactive_prefix, '/'];
        segmentations_filtered_location = [options.tmp_results_location, accumulated_directory_options_prefix, 'segmentations-filtered', segmentation_filtering_interactive_prefix, '/'];
        accumulated_directory_options_prefix = [accumulated_directory_options_prefix, 'segmentation-filtering', segmentation_filtering_interactive_prefix, '_'];
    end
else
    segmentations_filtered_location = segmentations_location;
end

% warning('Adding two-point-related code!!!')
% Add whether using two-point synapses to the prefix:
if options.use_two_point_synapses
    if should_use_short_directory_names
        accumulated_directory_options_prefix = [accumulated_directory_options_prefix, '2pt_'];
    else
        accumulated_directory_options_prefix = [accumulated_directory_options_prefix, 'two-point_'];
    end
end

% Add whether using landmark smoothing to the prefix:
if should_use_short_directory_names
    accumulated_directory_options_prefix = [accumulated_directory_options_prefix, sprintf('ldmksmmeth_%s_', options.landmark_temporal_smoothing_method)];
else
    accumulated_directory_options_prefix = [accumulated_directory_options_prefix, sprintf('landmark-smoothing-method_%s_', options.landmark_temporal_smoothing_method)];
end
switch options.landmark_temporal_smoothing_method
    case 'none'
        
    case 'lowess'
        % Smooth all frames' alignments:
        if should_use_short_directory_names
            accumulated_directory_options_prefix = [accumulated_directory_options_prefix, sprintf('smprpn_%.2f_', options.landmark_temporal_smoothing_proportion)];
        else
            accumulated_directory_options_prefix = [accumulated_directory_options_prefix, sprintf('smooth-proportion_%.2f_', options.landmark_temporal_smoothing_proportion)];
        end
    case 'lowess-for-outliers'
      if should_use_short_directory_names
        accumulated_directory_options_prefix = [accumulated_directory_options_prefix, sprintf('smprpn_%.2f_', options.landmark_temporal_smoothing_proportion)];
        % % Smooth only alignments of frames with this ratio of absolute intensity error to total intensity or greater:
        % Not clear yet how to implement this[t_cell_info] = tcell_automatic_detect_synapse(t_cell_info, options) filtering, probably should use training data instead of hand-tuning:
        accumulated_directory_options_prefix = [accumulated_directory_options_prefix, sprintf('outlthr_%.2f_', options.landmark_temporal_smoothing_outlier_relative_error)];
      else
        accumulated_directory_options_prefix = [accumulated_directory_options_prefix, sprintf('smooth-proportion_%.2f_', options.landmark_temporal_smoothing_proportion)];
        % % Smooth only alignments of frames with this ratio of absolute intensity error to total intensity or greater:
        % Not clear yet how to implement this filtering, probably should use training data instead of hand-tuning:
        accumulated_directory_options_prefix = [accumulated_directory_options_prefix, sprintf('outlier-threshold_%.2f_', options.landmark_temporal_smoothing_outlier_relative_error)];
      end
      error('Not yet implemented')
    otherwise
        error('Unrecognized option value landmark_temporal_smoothing_method = %s', options.landmark_temporal_smoothing_method)
end

if should_use_short_directory_names
    alignments_location = [options.tmp_results_location, accumulated_directory_options_prefix, 'algt/'];
    deformations_location = [options.tmp_results_location, accumulated_directory_options_prefix, 'dfmn/'];
else
    alignments_location = [options.tmp_results_location, accumulated_directory_options_prefix, 'alignments/'];
    deformations_location = [options.tmp_results_location, accumulated_directory_options_prefix, 'deformations/'];
end

if should_use_short_directory_names
    accumulated_directory_options_prefix = [accumulated_directory_options_prefix, sprintf('intymrphmeth_%s_', options.intensity_morphing_inside_template_method)];
else
    accumulated_directory_options_prefix = [accumulated_directory_options_prefix, sprintf('intensity-morphing-method_%s_', options.intensity_morphing_inside_template_method)];
end
switch options.intensity_morphing_inside_template_method
    case 'none'
    case 'lddmm'
        error('Not yet implemented')
    case 'ants'
        warning('Being implemented')
    otherwise
        error('Unrecognized option value intensity_morphing_inside_template_method = %s', options.intensity_morphing_inside_template_method)
end

if should_use_short_directory_names
    deformations_aligned_location = [options.tmp_results_location, accumulated_directory_options_prefix, 'dfmnalgd/'];
    models_location = [options.tmp_results_location, accumulated_directory_options_prefix, 'mdl/'];
else
    deformations_aligned_location = [options.tmp_results_location, accumulated_directory_options_prefix, 'deformations-aligned/'];
    models_location = [options.tmp_results_location, accumulated_directory_options_prefix, 'models/'];
end

illustration_location = [options.tmp_results_location, 'illustrations/'];


% xruan 07/18/2015
% add the paths for some figures
mkdir(options.results_location)
mkdir(regions_location)
mkdir(segmentations_location)
if options.use_segmentation_filtering
    mkdir(segmentation_evaluations_location)
end
mkdir(segmentations_filtered_location)
mkdir(alignments_location)
mkdir(deformations_location)
mkdir(deformations_aligned_location)
mkdir(models_location)
mkdir(cell_positions_location)
mkdir(illustration_location)
mkdir(temp_location)

t_cell_info.path_info = struct();
t_cell_info.path_info.regions_location = regions_location;
t_cell_info.path_info.segmentations_location = segmentations_location;
if options.use_segmentation_filtering
    t_cell_info.path_info.segmentation_evaluations_location = segmentation_evaluations_location;
end
t_cell_info.path_info.segmentations_filtered_location = segmentations_filtered_location;
t_cell_info.path_info.alignments_location = alignments_location;
t_cell_info.path_info.deformations_location = deformations_location;
t_cell_info.path_info.deformations_aligned_location = deformations_aligned_location;
t_cell_info.path_info.models_location = models_location;
t_cell_info.path_info.cell_positions_location = cell_positions_location;
t_cell_info.path_info.illustration_location = illustration_location;
t_cell_info.path_info.temp_location = temp_location;

if options.use_WAVE2_manual_seg
 t_cell_info.path_info.WAVE2_active_root_dir  = options.WAVE2_active_root_dir;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure out how large the window surrounding a synapse point should be and create a template cell shape of similar size

[t_cell_info] = tcell_get_template_info(t_cell_info);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Colormap with which to display sensor intensities:

t_cell_info.preprocessing_info.image_channels = options.image_channels;
% disable setting of colormap in the main pipeline 
% illustration_colormap_to_use_name = 'pmkmp_cubic';
% switch illustration_colormap_to_use_name   
%   case 'pmkmp_cubic'
%     illustration_colormap_to_use = pmkmp(1024, 'CubicL'); 
% end 

t_cell_info.plotting_info = struct();

% warning('IDK if this works yet!')
env_display = getenv('DISPLAY');
env_display = regexprep(env_display, '[\r\n]', '');
if isempty(env_display)
env_display = '';
end
% t_cell_info.plotting_info.running_interactively = ~isempty(env_display);
% t_cell_info.plotting_info.illustration_colormap_to_use_name = illustration_colormap_to_use_name;
% t_cell_info.plotting_info.illustration_colormap_to_use = illustration_colormap_to_use;
% t_cell_info.plotting_info.grayscale_color_mapping_function = @grayscale_color_mapping_function;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set random seed: 

s = RandStream('mt19937ar', 'Seed', 892734);
% xruan 06/2015 
% update the code, setDefaultStream is out of date.
RandStream.setGlobalStream(s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Infer synapse based on only annotation of one point
if isfield(options, 'infer_synapses') && options.infer_synapses
  inferred_synapse_location = [options.temporary_results, '/synapse_inference/'];
  mkdir(inferred_synapse_location);
  t_cell_info.path_info.inferred_synapse_location = inferred_synapse_location;
  [t_cell_info] = tcell_automatic_detect_synapse(t_cell_info, options);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Produce lists of synapse files, images, and individual synapses, frames, and cells:

% [t_cell_info] = tcell_get_synapse_info(t_cell_info);
[t_cell_info] = tcell_get_synapse_info_new(t_cell_info);

% set up default model_prefix if it is empty
if isempty(options.model_prefix)
  run_relative_path_list = t_cell_info.synapse_info.run_relative_path_list;      
  options.model_prefix = [run_relative_path_list{1}(1 : end - 4), '_'];
  t_cell_info.options.model_prefix = options.model_prefix;
end

% get segmentation parameters 
[t_cell_info] = tcell_get_segmentation_info(t_cell_info);
% get parameters for morphing using diffeomorphic model
[t_cell_info] = tcell_get_morphing_info(t_cell_info);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options.t_cell_info = t_cell_info;
cell_images = t_cell_info.synapse_info.image_name_run_list;

end
