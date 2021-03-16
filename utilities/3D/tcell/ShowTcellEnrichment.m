function [enrichment_means, enrichment_stds, all_enrichments, all_timepoints, h] = ShowTcellEnrichment(model_filenames, param)
% The function shows how much of the protein is enriched in the
% (approximate) synapse region at various time points.
% The input is the model_filenames. Currently only support single protein
% enrichment plot.
%
%
% Author: Xiongtao Ruan
% Date: Sep. 2016, Jan. 2017
% 11/05/2018 add param to control the behavior of the function, and add new
% functionality for using user-defined enrichment region, and enrichment of
% one region over another region.
% 11/6/2018 add outputs: mean_e
% 11/08/2018 add enrichment region type
% 11/14/2018 add enrichment region over 90% fluorescence region
% 11/26/2018 add option for errorbar as sd or sem; choose to save result if
% providing a filename
% 01/29/2019 fix bug for single cell enrichment analysis, and also enable
% spaces in the path to the models
% 03/01/2021 R.F.Murphy change interpolation from CubicL to PCHIP due to
% Matlab deprecating CubicL
% 03/02/2021 R.F.Murphy undo change to CubicL and add documentation for options

% values for options (all prefixed by model.tcell.)
%
% error_bar_type - 'sd' or 'sem' - default 'sd'
% save_result_filename - if specified, save enrichment values to .csv file - default none
% enrichment_region_percentiles - default [90]
% should_use_global_enrichment_region - true or false - 
%   define region using enrichment_region_percentiles - default true
% should_use_user_defined_enrichment_region' - true or false - default false
% enrichment_over_certain_region - true or false - 
%   define both a numerator region and a denominator region - default false
% exclude_enrichment_region - true or false - 
%   don't include enrich regionin denominator of percentages - default false
% enrichment_region_type - 'cylinder' or 'ring' or 'top_fluorescence' - default 'cylinder'
% enrichment_bottom_region_type' - 'cylinder' or 'ring' - default 'cylinder'
% enrichment_bottom_region_percentiles - , 90, ...
% region_radius - numeric - default 8
% region_thickness - numeric - default 4
% region_start_ind - numeric - default 1
% region_2_radius - numeric - for defining ring region, default 8 / 3
% region_2_thickness - numeric - for defining ring region - default 4
% region_2_start_ind - numeric - for defining ring region - default 1

debug = true;

if nargin < 1
    error('Models argument is missing');
end

% create the parameter structure and include default structure.
if nargin < 2
    param = struct();
end

% Filter the options to only the one set for tcell model
if isfield(param, 'model') && isfield(param.model, 'tcell')
  param = param.model.tcell;
else
  param = struct();
end

param = process_options_structure(struct('error_bar_type', 'sd', ...
                                         'save_result_filename', '', ...
                                         'enrichment_region_percentiles', [90], ...
                                         'should_use_per_cell_max_enrichment', false, ...
                                         'should_use_per_cell_enrichment_region', false, ...
                                         'should_use_actin_full_stimulus_enrichment_region', false, ...
                                         'should_use_global_enrichment_region', true, ...
                                         'should_use_user_defined_enrichment_region', false, ...
                                         'enrichment_over_certain_region', false, ...
                                         'exclude_enrichment_region', false, ...
                                         'enrichment_region_type', 'cylinder', ...
                                         'enrichment_bottom_region_type', 'cylinder', ...
                                         'enrichment_bottom_region_percentiles', 90, ...
                                         'region_radius', 8, ...
                                         'region_thickness', 4, ...
                                         'region_start_ind', 1, ...
                                         'region_2_radius', 8 / 3, ...
                                         'region_2_thickness', 4, ...
                                         'region_2_start_ind', 1), param);


enrichment_region_percentiles = param.enrichment_region_percentiles;
should_use_per_cell_max_enrichment = param.should_use_per_cell_max_enrichment;
should_use_per_cell_enrichment_region = param.should_use_per_cell_enrichment_region;
should_use_actin_full_stimulus_enrichment_region = param.should_use_actin_full_stimulus_enrichment_region;
should_use_global_enrichment_region = param.should_use_global_enrichment_region;
should_use_user_defined_enrichment_region = param.should_use_user_defined_enrichment_region;
enrichment_over_certain_region = param.enrichment_over_certain_region;
exclude_enrichment_region = param.exclude_enrichment_region;
enrichment_region_type = param.enrichment_region_type;
enrichment_bottom_region_type = param.enrichment_bottom_region_type;
enrichment_bottom_region_percentiles = param.enrichment_bottom_region_percentiles;


% errorbar type be 'sd' or 'sem'
error_bar_type = param.error_bar_type;
% if save_result_filename not empty, save the enrichment values for all
% time points
save_result_filename = param.save_result_filename;

if should_use_user_defined_enrichment_region
    region_radius = param.region_radius;
    region_thickness = param.region_thickness;
    region_start_ind = param.region_start_ind;
end

if enrichment_over_certain_region
    region_2_radius = param.region_2_radius;
    region_2_thickness = param.region_2_thickness;
    region_2_start_ind = param.region_2_start_ind;
end

% 01/29/2019 here we assume model_filenames is just a single string that
% may contains space
if ~iscell(model_filenames)
  if contains(model_filenames, ' ')
      model_filenames_temp = strrep(model_filenames, '\ ', ' ');
      model_filenames = strrep(model_filenames_temp, ' ', '\ ');
  end
  model_file_list = ml_ls(model_filenames);
else
    model_file_list = model_filenames;
end
model_numbers = numel(model_file_list);

try
    load(model_file_list{1});
catch
    error('Unable to load model file!');
end

t_cell_info = model.proteinModel.t_cell_info;

% Variables referenced more than once can be copied to local variables:
master_script_options = t_cell_info.options;
regions_location = t_cell_info.path_info.regions_location;
segmentations_filtered_location = t_cell_info.path_info.segmentations_filtered_location;
alignments_location = t_cell_info.path_info.alignments_location;
deformations_location = t_cell_info.path_info.deformations_location;
models_location = t_cell_info.path_info.models_location;

% illustration_colormap_to_use = t_cell_info.plotting_info.illustration_colormap_to_use;

illustration_colormap_to_use_name = 'pmkmp_cubic';
switch illustration_colormap_to_use_name
    case 'pmkmp_cubic'
        illustration_colormap_to_use = pmkmp(1024, 'CubicL');
end

condition_name_abbreviation_function = t_cell_info.abbreviation_info.condition_name_abbreviation_function;
model_types = t_cell_info.model_type_info.model_types;
number_model_types = t_cell_info.model_type_info.number_model_types;
model_representation_functions = t_cell_info.model_type_info.model_representation_functions;
model_reconstruction_functions = t_cell_info.model_type_info.model_reconstruction_functions;

model_reconstruction_function = model_reconstruction_functions{1};
model_type = model_types{1};

template_centroid = t_cell_info.template_info.template_centroid;
template_synapse = t_cell_info.template_info.template_synapse;
template_image = t_cell_info.template_info.template_image >= 0.5;
template_crop_function = t_cell_info.template_info.template_crop_function;

window_size_2d = t_cell_info.preprocessing_info.window_size_2d;
window_center_2d = t_cell_info.preprocessing_info.window_center_2d;
cropped_size = t_cell_info.preprocessing_info.cropped_size;

segmentation_rasterization_function = t_cell_info.segmentation_info.segmentation_rasterization_function;

should_compute_individual_frame_illustrations = false;

should_skip_completed_images = false;
current_mean_intensity_limit_excluded_size = 0;
intensity_limit_function = @(given_image, given_limits)given_image .* (given_image < given_limits(2)) .* (given_image > given_limits(1)) + given_limits(1) .* (given_image <= given_limits(1)) + given_limits(2) .* (given_image >= given_limits(2));

model_relative_times = zeros(model_numbers, 1);
model_mean = [];
model_count = 0;

% compute model mean for all time points
for i = 1 : model_numbers
    cur_model_filename = model_file_list{i};
    try
        load(cur_model_filename);
    catch
        error('Unable to load model file!');
    end
    current_data = model.proteinModel.current_data;
    model_relative_times(i) = model.proteinModel.relative_time;

    if isempty(model_mean)
        model_mean = zeros(size(current_data(1, :)));
    end
     model_mean = (model_mean * model_count + sum(current_data, 1)) / (model_count + size(current_data, 1));
     model_count = model_count + size(current_data, 1);
end

% sort model name so that they are from small to large time points
[model_relative_times, inds] = sort(model_relative_times);
model_file_list = model_file_list(inds);

% warning('Why doesn''t reconstruction_function get modified before this call???'), % model_mean_reconstructed = reconstruction_function(model_mean);
model_mean_reconstructed = model_reconstruction_function(model_mean);
% Take the top 2% of voxels from the overall mean image as the region of interest for any enrichment (looking at the data, assumes that all interesting enrichment is at the synapse and ):
% enrichment_region_percentile = 98;
% model_mean_enrichment_region = model_mean_reconstructed >= prctile(model_mean_reconstructed(template_image), enrichment_region_percentile);
model_mean_enrichment_region_function = @(given_percentile)model_mean_reconstructed >= prctile(model_mean_reconstructed(template_image), given_percentile);


% if use user defined region, calculate the region
if should_use_user_defined_enrichment_region
    if enrichment_over_certain_region
        cur_enrichment_region_type = enrichment_bottom_region_type;
    else
        cur_enrichment_region_type = enrichment_region_type;
    end

    switch cur_enrichment_region_type
        case 'cylinder'
            user_defined_region = Tcell_user_defined_synapse_region(template_image, region_radius, region_thickness, region_start_ind);
        case 'ring'
            user_defined_region_1_1 = Tcell_user_defined_synapse_region(template_image, region_radius(1), region_thickness(1), region_start_ind(1));
            user_defined_region_1_2 = Tcell_user_defined_synapse_region(template_image, region_radius(2), region_thickness(2), region_start_ind(2));
            user_defined_region = user_defined_region_1_1 - user_defined_region_1_2 > 0;
        case 'top_fluorescence'
            user_defined_region = model_mean_enrichment_region_function(enrichment_bottom_region_percentiles);
    end
end

if enrichment_over_certain_region
    switch enrichment_region_type
        case 'cylinder'
            user_defined_region_2 = Tcell_user_defined_synapse_region(template_image, region_2_radius, region_2_thickness, region_2_start_ind);
        case 'ring'
            user_defined_region_2_1 = Tcell_user_defined_synapse_region(template_image, region_2_radius(1), region_2_thickness(1), region_2_start_ind(1));
            user_defined_region_2_2 = Tcell_user_defined_synapse_region(template_image, region_2_radius(2), region_2_thickness(2), region_2_start_ind(2));
            user_defined_region_2 = user_defined_region_2_1 - user_defined_region_2_2 > 0;
    end
end

% 11/05/2018 set for different conditions.
if should_use_user_defined_enrichment_region
    if enrichment_over_certain_region
        model_mean_enrichment_regions = {user_defined_region_2};
    else
        model_mean_enrichment_regions = {user_defined_region};
    end
else
    model_mean_enrichment_regions = arrayfun(model_mean_enrichment_region_function, enrichment_region_percentiles, 'UniformOutput', false);
end

if false
  imshow(reshape_contrast([model_mean_reconstructed; model_mean_reconstructed .* model_enrichment_region], -1))
  % pause
  whos model_mean model_mean_reconstructed, model_count
  keyboard
end

% Write the overall mean image and the above-threshold portion as an image:
for enrichment_region_percentile_index = 1:length(model_mean_enrichment_regions)
  enrichment_region_percentile = enrichment_region_percentiles(enrichment_region_percentile_index);
    model_mean_enrichment_region = model_mean_enrichment_regions{enrichment_region_percentile_index};
  image_to_show = [];
  % image_to_show = [contrast_stretch(model_mean_reconstructed); model_mean_enrichment_region];
  % image_to_show = [contrast_stretch(model_mean_reconstructed); contrast_stretch(model_mean_enrichment_region * 3 + template_image)];
  % image_to_show = [image_to_show; contrast_stretch(model_mean_reconstructed)];
  image_to_show = [image_to_show; contrast_stretch(template_crop_function(model_mean_reconstructed))];
  % image_to_show = [image_to_show; model_mean_enrichment_region_function(enrichment_region_percentile) * .675 + template_image * .2];
  % image_to_show = [image_to_show; model_mean_enrichment_region_function(enrichment_region_percentile) * .675 + template_image * .2];
  % image_to_show = [image_to_show; model_mean_enrichment_region * .675 + template_image * .2];
  % image_to_show = [image_to_show; template_crop_function(model_mean_enrichment_region * .675 + template_image * .2)];

  if should_use_user_defined_enrichment_region
    if enrichment_over_certain_region
      % image_to_show = [image_to_show; template_crop_function(user_defined_region .* 0.675 + model_mean_enrichment_region * .375 + template_image * .2)];
      image_to_show = [image_to_show; template_crop_function(model_mean_enrichment_region * .675 + template_image * .2)];
      image_to_show = [image_to_show; template_crop_function(user_defined_region .* 0.675 + template_image * .2)];
    else
      image_to_show = [image_to_show; template_crop_function(model_mean_enrichment_region * .675 + template_image * .2)];
    end
  else
    image_to_show = [image_to_show; template_crop_function(model_mean_enrichment_region * .675 + template_image * .2)];
  end

  % image_to_show = reshape_2d(image_to_show, -1);
  image_to_show = reshape_2d(image_to_show, 1);
  % Add color bar:
  image_to_show = [image_to_show, repmat(linspace(1, 0, size(image_to_show, 1)).', [1, size(template_crop_function(template_image), 2)])];
  % Recolor the image:
  image_to_show = ind2rgb(floor(contrast_stretch(image_to_show) .* (size(illustration_colormap_to_use, 1) - 1)) + 1, illustration_colormap_to_use);
  % Write:
  figure,
  imshow(image_to_show, []);
  image_filename = [sprintf('%smodel_mean_enrichment_region_mdl-type_%s_enr-prct%d', master_script_options.model_prefix, model_type, enrichment_region_percentile)];
  image_full_filename = [ image_filename, '.png'];
  imwrite(image_to_show, image_full_filename)
  fprintf('Wrote figure %s.png\n', image_filename)
end

% enrichment for each time point.

enrichment_types = {};
enrichment_type_computation_functions = {};
enrichment_type_mean_computation_functions = {};

% 11/05/2018 add enrichment by user defined region.
if should_use_user_defined_enrichment_region
    if enrichment_over_certain_region
      enrichment_types{end + 1} = 'user\_defined\_two\_regions';
      user_defined_region = user_defined_region >= 0.5;
      user_defined_region_2 = user_defined_region_2 >= 0.5;
      if exclude_enrichment_region
        enrichment_type_computation_functions{end + 1} = @(x) mean(x(user_defined_region_2)) ./ mean(x((user_defined_region - user_defined_region_2) >= .5));
      else
        enrichment_type_computation_functions{end + 1} = @(x) mean(x(user_defined_region_2)) ./ mean(x(user_defined_region));
      end
    else
      enrichment_types{end + 1} = 'user\_defined\_region';
      user_defined_region = user_defined_region >= 0.5;
      if exclude_enrichment_region
        enrichment_type_computation_functions{end + 1} = @(x) mean(x(user_defined_region)) ./ mean(x((template_image - user_defined_region) >= 0.5));
      else
        enrichment_type_computation_functions{end + 1} = @(x) mean(x(user_defined_region)) ./ mean(x(template_image));
      end
    end
else
    if should_use_per_cell_max_enrichment
      % Max / mean intensity:
      enrichment_types{end + 1} = 'max';
      enrichment_type_computation_functions{end + 1} = @(x)max(x(template_image >= .5)) ./ mean(x(template_image >= .5));
    end
    % for enrichment_region_percentile = enrichment_region_percentiles
    for enrichment_region_percentile_index = 1:length(enrichment_region_percentiles)
      enrichment_region_percentile = enrichment_region_percentiles(enrichment_region_percentile_index);
      model_mean_enrichment_region = model_mean_enrichment_regions{enrichment_region_percentile_index};
      if should_use_per_cell_enrichment_region
        % Mean intensity in per-cell enrichment region / mean intensity:
        enrichment_types{end + 1} = sprintf('celltop%03d', enrichment_region_percentile);
        enrichment_type_computation_functions{end + 1} = @(x)mean(x(x >= prctile(x(template_image >= .5), enrichment_region_percentile))) ./ mean(x(template_image >= .5));
      end
      if should_use_global_enrichment_region
        % Mean intensity in mean's enrichment region / mean intensity:
        enrichment_types{end + 1} = sprintf('meantop%03d', enrichment_region_percentile);
        % enrichment_type_computation_functions{end + 1} = @(x)mean(x(model_mean_enrichment_region_function(enrichment_region_percentile))) ./ mean(x(template_image >= .5));
          if exclude_enrichment_region
            enrichment_type_computation_functions{end + 1} = @(x)mean(x(model_mean_enrichment_region)) ./ mean(x((template_image - model_mean_enrichment_region) >= .5));
          else
            enrichment_type_computation_functions{end + 1} = @(x)mean(x(model_mean_enrichment_region)) ./ mean(x(template_image >= .5));
          end
      end
      % xruan 12/11/2015
      if should_use_actin_full_stimulus_enrichment_region
          actin_model_mean_enrichment_region = actin_model_mean_enrichment_regions{enrichment_region_percentile_index};
          enrichment_types{end + 1} = sprintf('ActinFullMeanTop%03d', enrichment_region_percentile);
          enrichment_type_computation_functions{end + 1} = @(x)mean(x(actin_model_mean_enrichment_region)) ./ mean(x(template_image >= .5));
      end
    end
end

number_enrichment_types = length(enrichment_types);
all_enrichments_mean = cell(number_enrichment_types);
all_enrichments_std = cell(number_enrichment_types);

for enrichment_type_index = 1:number_enrichment_types
    current_enrichment_type = enrichment_types{enrichment_type_index};
    current_enrichment_type_computation_function = enrichment_type_computation_functions{enrichment_type_index};

    current_enrichments = [];
    current_all_enrichments = {};
    current_enrichment_means = [];
    current_enrichment_stds = [];

    for model_index = 1 : model_numbers
        cur_model_filename = model_file_list{model_index};
        try
            load(cur_model_filename);
        catch
            error('Unable to load model file!');
        end
        current_data = model.proteinModel.current_data;
        current_enrichments = zeros(size(current_data, 1), 1);
        for datum_index = 1 : size(current_data, 1)
          datum_reconstruction = model_reconstruction_function(current_data(datum_index, :));
          % current_enrichments(end + 1) = max(datum_reconstruction(template_image >= .5)) ./ mean(datum_reconstruction(template_image >= .5));
          datum_enrichment = current_enrichment_type_computation_function(datum_reconstruction);
          current_enrichments(datum_index) = datum_enrichment;
        end
        current_all_enrichments{model_index} = current_enrichments;
        current_enrichment_mean = mean(current_enrichments);
        current_enrichment_std = std(current_enrichments);
        switch error_bar_type
            case 'sd'
            case 'sem'
                current_enrichment_std = current_enrichment_std / sqrt(numel(current_enrichments));
        end
        current_enrichment_means(model_index) = current_enrichment_mean;
        current_enrichment_stds(model_index) = current_enrichment_std;
    end
    all_enrichments_mean{enrichment_type_index} = current_enrichment_means;
    all_enrichments_std{enrichment_type_index} = current_enrichment_stds;
    figure,
    % by default we only calculate one type of enrichment
    h = plot(model_relative_times, current_enrichment_means);
    errorbar(model_relative_times, current_enrichment_means, current_enrichment_stds);
    title(sprintf('Enrichment plot type : %s', current_enrichment_type));
    ylabel('Enrichment');
    xlabel('relative times');
    figure_filename = [sprintf('%smodel_mean_enrichment_plot_mdl-type_%s_enr-prct%d', master_script_options.model_prefix, model_type, enrichment_region_percentile)];
    figure_full_filename = [ figure_filename, '.png'];
    saveas(gcf, figure_full_filename)
    fprintf('Wrote figure %s.png\n', figure_filename)
end

% outputs
enrichment_means = current_enrichment_means;
enrichment_stds = current_enrichment_stds;
all_enrichments = current_all_enrichments;
all_timepoints = model_relative_times;

% if the filename string is not empty save the results.
if ~isempty(save_result_filename)
    if numel(save_result_filename) < 4 || ~strcmp(save_result_filename(end -3 : end), '.csv')
        save_result_filename = [save_result_filename, '.csv'];
    end
    fprintf('Write enrichment result to %s ...', save_result_filename);
    max_length = max(cellfun(@numel, all_enrichments));
    all_enrichments_cell = cell(max_length, numel(all_timepoints));

    for i = 1 : numel(all_timepoints)
        cur_enrichments = all_enrichments{i};
        all_enrichments_cell(1 : numel(cur_enrichments), i) = mat2cell(cur_enrichments, ones(numel(cur_enrichments), 1));
        all_enrichments_cell(numel(cur_enrichments) + 1 : end, i) = repmat({''}, max_length - numel(cur_enrichments), 1);
    end
    variable_names = arrayfun(@(x) sprintf('reltime %d', x), all_timepoints, 'uniformoutput', false);
    % T = cell2table(all_enrichments_cell, 'VariableNames', variable_names)

    fid = fopen(save_result_filename, 'w');
    fid = fopen(save_result_filename, 'wt');
    header_string = strjoin(variable_names, ',');
    fwrite(fid, header_string);
    fprintf(fid, '\n');

    for i = 1 : max_length
        cur_line = all_enrichments_cell(i, :);
        for k = 1 : numel(cur_line)
            if k > 1
                fprintf(fid, ',');
            end
            if ischar(cur_line{k})
                fwrite(fid, cur_line{k});
            else
                fprintf(fid, '%d', cur_line{k});
            end
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
end


end
