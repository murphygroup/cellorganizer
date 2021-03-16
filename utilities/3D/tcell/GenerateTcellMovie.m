function GenerateTcellMovie(models_1, models_2, models_3, protein_names)
% The script is for generating a movie showing protein intensity change
% across different time points. 
% The first three input are model names, each is a set of models for a
% protein. And if there are multiple proteins, their time points need to be
% exactly match each other. 
% The fourth input is optional, it define the protein names use a cell
% array. If there is no such input, the function will use the model prefix,
% that is the string before '_reltime_*mat' of the first model name for the
% protein. 
% The default video dir is ./video relative to current directory. The video
% is in .avi format. And the intermediate images used for generated for the
% video is also there. 
% By default, we show the top 25 intensity of the frame. For advanced user,
% Other threshold or visualization methods could be used. 
%
% Author: Xiongtao Ruan
% Date: Jan. 2017
%
%


debug = ~true;
if debug
    models_1 = './models/LAT_reltime_*.mat';
    models_2 = './models/Ezrin_reltime_*.mat';   
    % models_2 = [];
    models_3 = [];
    protein_names = {'LAT', 'Ezrin'};
    % protein_names = {'LAT'};
end
if ~debug
    if nargin < 1
        error('Please provide at least one set of models!');
    end
    if nargin < 2
         models_2 = [];
    end
    if nargin < 3
         models_3 = [];
    end
    if nargin < 4
        protein_names = {};
    end
end

models_list = {models_1, models_2, models_3};
models_list(cellfun(@isempty, models_list)) = [];
model_relative_times_mat = [];
for i = 1 : numel(models_list)
    cur_model_list = ml_ls(models_list{i});
    [model_relative_times, model_name_list] = extract_model_relative_times(cur_model_list);
    models_list{i} = model_name_list;
    model_relative_times_mat = [model_relative_times_mat, model_relative_times];
end

if isempty(protein_names)
    for i = 1 : numel(models_list)
        current_first_model_name = models_list{i}{1};
        [~, current_first_model_name] = fileparts(current_first_model_name);
        current_protein_name = regexp(current_first_model_name, '(.*)_reltime', 'tokens');
        if ~isempty(current_protein_name)
            protein_names{i} = current_protein_name{1}{1};
        end
    end
end

relative_time_diff = pdist(model_relative_times_mat');
if any(relative_time_diff(:) > 0) 
    error('All proteins need to have the same set of relative time points.');
end
all_relative_times = model_relative_times_mat(:, 1);

num_models = numel(models_list);

try 
    load(models_list{1}{1});
catch
    error('Unable to load model file!');
end

t_cell_info = model.proteinModel.t_cell_info;

timepoints_to_include = all_relative_times;

% Variables referenced more than once can be copied to local variables:
master_script_options = t_cell_info.options;
regions_location = t_cell_info.path_info.regions_location;
segmentations_filtered_location = t_cell_info.path_info.segmentations_filtered_location;
alignments_location = t_cell_info.path_info.alignments_location;
deformations_location = t_cell_info.path_info.deformations_location;
models_location = t_cell_info.path_info.models_location;
%   figures_location = t_cell_info.path_info.figures_location;
%   videos_location = t_cell_info.path_info.videos_location;
master_script_options.video_prefix = strjoin(protein_names, '_');

% removed the option before. 
% illustration_colormap_to_use = t_cell_info.plotting_info.illustration_colormap_to_use;
illustration_colormap_to_use_name = 'pmkmp_cubic';
switch illustration_colormap_to_use_name
    case 'pmkmp_cubic'
        illustration_colormap_to_use = pmkmp(1024, 'CubicL'); 
end 

% xruan 08/17/2015
% add the usage of condition name abbreviation

model_types = t_cell_info.model_type_info.model_types;
number_model_types = t_cell_info.model_type_info.number_model_types;
model_representation_functions = t_cell_info.model_type_info.model_representation_functions;
model_reconstruction_functions = t_cell_info.model_type_info.model_reconstruction_functions;

template_centroid = t_cell_info.template_info.template_centroid;
template_synapse = t_cell_info.template_info.template_synapse;
template_image = t_cell_info.template_info.template_image;
template_crop_function = t_cell_info.template_info.template_crop_function;

window_size_2d = t_cell_info.preprocessing_info.window_size_2d;
window_center_2d = t_cell_info.preprocessing_info.window_center_2d;
cropped_size = t_cell_info.preprocessing_info.cropped_size;

segmentation_rasterization_function = t_cell_info.segmentation_info.segmentation_rasterization_function;

current_mean_intensity_limit_excluded_size = 0;
% current_mean_intensity_limit_excluded_size = 1e-3;
% current_mean_intensity_limit_excluded_size = 5e-3;
intensity_limit_function = @(given_image, given_limits)given_image .* (given_image < given_limits(2)) .* (given_image > given_limits(1)) + given_limits(1) .* (given_image <= given_limits(1)) + given_limits(2) .* (given_image >= given_limits(2));

if master_script_options.illustration_uses_histogram_equalization
    % histogram_equalization_function = @(given_image)masked_histogram_equalization(given_image, template_image > .5);
    histogram_equalization_function = @(given_image, given_limits)masked_histogram_equalization(given_image, (given_image > given_limits(1)) & (given_image < given_limits(2)));
else
    histogram_equalization_function = @(given_image, given_limits)given_image;
end

if master_script_options.illustration_uses_logarithmic_scale
    intensity_warping_function = @(given_image, given_limits)log(given_image);
else
    intensity_warping_function = @(given_image, given_limits)given_image;
end
  

for model_type_index = 1:number_model_types
    model_type = model_types{model_type_index};
    representation_function = model_representation_functions{model_type_index};
    reconstruction_function = model_reconstruction_functions{model_type_index};

    model_type

    % Combine some sensors into multi-channel videos:
    % Each entry in this cell array is a 3x2 cell array containing condition as the first column and sensor as the second (condition is assumed constant across channels, and 3 channels assumed, for the code below):
    % video_channels = {...
      % {'Full Stimulus', 'ARP3'; 'Full Stimulus', 'CPalpha1'; 'Full Stimulus', 'WASP'}; ...
      % {'Full Stimulus', 'Cofilin'; 'Full Stimulus', 'MRLC'; 'Full Stimulus', 'WASP'}; ...
      % };
    video_channels = {};
    video_titles = {};

    should_videos_be_colorblind_friendly = false;
    % should_videos_be_colorblind_friendly = true;

    % Then go through each set of two (should_videos_be_colorblind_friendly) or three clusters:
    current_channel_index = 1;
    current_channels = cell(1, 3);

    % video_channels
    % video_titles
    % beep, keyboard
    
    video_channels{1} = models_list';


    videos_location = [pwd, '/' 'video/'];

    if ~exist(videos_location, 'dir')
        mkdir(videos_location);
    end

    % current_video_frame_index = 1;
    for video_index = 1:length(video_channels)
        current_video_channels = video_channels{video_index};
        current_video_title = '';
        num_timepoints_to_include = numel(timepoints_to_include);
        for all_relative_time_index = 1:num_timepoints_to_include
            relative_time = timepoints_to_include(all_relative_time_index);

            relative_time_image_rgb = [];
            current_video_frame_index = (video_index - 1) * num_timepoints_to_include + (all_relative_time_index - 1) + 1;

            for channel_index = 1:size(current_video_channels, 1)
                current_model_names = current_video_channels{channel_index, 1};

                % xruan 08/17/2015 
                % use condition abbreviation name as did in
                % master_script_generate_modeififl_figures
                model_index = all_relative_time_index;
                current_model_name = current_model_names{model_index};
                image_filename = regexp(current_model_name, '([^/]+).mat', 'tokens');
                image_filename = image_filename{1}{1};

                % image_filename = sprintf('%s_%s_reltime% 3d_%s', current_condition, current_sensor, relative_time, model_type);
                % image_filename = sprintf('%s%s_%s_frame%03d_reltime% 3d_%s', master_script_options.model_prefix, current_condition, current_sensor, all_relative_time_index, relative_time, model_type);

                grayscale_image_full_filename = [videos_location, image_filename, '_grayscale.png'];
                try 
                  a = load(current_model_name);
                catch
                  error('unable to load model')
                end
                current_mean = a.model.proteinModel.current_mean;
                current_mean_reconstructed = reconstruction_function(current_mean);

                % Produce the image for this condition-sensor combination:

                current_mean_intensity_limits = [0, quantile(current_mean_reconstructed(template_image > .5), 1 - current_mean_intensity_limit_excluded_size)];


                current_mean_cropped = template_crop_function(current_mean_reconstructed);
                current_mean_cropped = intensity_limit_function(current_mean_cropped, current_mean_intensity_limits);
                current_mean_cropped = histogram_equalization_function(current_mean_cropped, current_mean_intensity_limits);
                current_mean_cropped = intensity_warping_function(current_mean_cropped, current_mean_intensity_limits);

                % current_mean_cropped_zy = permute(current_mean_cropped, [2, 3, 1]);
                current_mean_cropped_zy = permute(current_mean_cropped, [3, 2, 1]);

                image_to_show = [reshape_contrast(current_mean_cropped, -1); reshape_contrast(current_mean_cropped_zy, -1)];
                % grayscale_image_full_filename = [figures_location, image_filename, '_grayscale.png'];
                grayscale_image_full_filename = [videos_location, image_filename, '_grayscale.png'];
                imwrite(image_to_show, grayscale_image_full_filename)

                if ~exist(grayscale_image_full_filename, 'file')
                    % beep, keyboard
                    continue
                end
                grayscale_image = double(imread(grayscale_image_full_filename));
                grayscale_image = contrast_stretch(grayscale_image);
                % grayscale_image = imresize(grayscale_image, master_script_options.videos_scale, 'bilinear');
                grayscale_image = imresize(grayscale_image, master_script_options.videos_scale, 'nearest');

                master_script_options.videos_contrast_enhancement_method = 'top_intensity25';
                current_video_title = '';

              switch master_script_options.videos_contrast_enhancement_method
                case 'exp-scale'
                  grayscale_image = exp(grayscale_image);
                  grayscale_image(~isfinite(grayscale_image)) = min(grayscale_image(isfinite(grayscale_image)));
                  grayscale_image = contrast_stretch(grayscale_image);

                case 'log-scale'
                  grayscale_image = log(grayscale_image);
                  grayscale_image(~isfinite(grayscale_image)) = min(grayscale_image(isfinite(grayscale_image)));
                  grayscale_image = contrast_stretch(grayscale_image);

                case 'threshold'
                  if ~false
                    % Debug info:
                    close all
                    figure,
                    imshow([grayscale_image; grayscale_image >= graythresh(grayscale_image(grayscale_image > min(grayscale_image(:))))], [])
                     beep
                    pause(1)
                    % keyboard
                  end
                  % grayscale_image = grayscale_image >= graythresh(grayscale_image(template_image >= .5));
                  % Hacky, should use a properly reshaped template_crop_function of template_image:        
                  grayscale_image = grayscale_image >= graythresh(grayscale_image(grayscale_image > min(grayscale_image(:))));
                  grayscale_image = double(grayscale_image);

               case 'top_intensity25'
                 percetile_val = 0.75;
                 grayscale_image(145 : 192, 231 : 276) = 0;
                 intensity_val = grayscale_image(grayscale_image > min(grayscale_image(:)));
                 intensity_val_sorted = sort(intensity_val);
                 intensity_cdf = cumsum(intensity_val_sorted) / sum(intensity_val_sorted);
                 [~, ind] = min(abs(intensity_cdf - percetile_val));
                 threshold = intensity_val_sorted(ind);
                 if false
                    % Debug info:
                    graythresh(grayscale_image(grayscale_image > min(grayscale_image(:))));
                    threshold
                    close all
                    figure,
                    imshow([grayscale_image; double(grayscale_image > threshold)], []);
                    beep
                    pause(1)
                    % keyboard
                 end

                 grayscale_image = double(grayscale_image > threshold);          
                case 'top_intensity10'
                 percetile_val = 0.9;
                 grayscale_image(145 : 192, 231 : 276) = 0;
                 intensity_val = grayscale_image(grayscale_image > min(grayscale_image(:)));
                 intensity_val_sorted = sort(intensity_val);
                 intensity_cdf = cumsum(intensity_val_sorted) / sum(intensity_val_sorted);
                 [~, ind] = min(abs(intensity_cdf - percetile_val));
                 threshold = intensity_val_sorted(ind);
                 if false
                    % Debug info:
                    graythresh(grayscale_image(grayscale_image > min(grayscale_image(:))));
                    threshold
                    close all
                    figure,
                    imshow([grayscale_image; double(grayscale_image > threshold)], []);
                    beep
                    pause(1)
                    % keyboard
                  end
                  grayscale_image = double(grayscale_image > threshold);

                case 'none'

                % case ''

                otherwise
                  error
              end

              if isempty(relative_time_image_rgb)
                relative_time_image_rgb = zeros([size(grayscale_image), 3 * ones(1, 3 - ndims(grayscale_image))]);
              end
              if should_videos_be_colorblind_friendly
                relative_time_image_rgb(:, :, channel_index) = grayscale_image(:, :, 1);
                if channel_index == 1
                  relative_time_image_rgb(:, :, 3) = grayscale_image(:, :, 1);
                end
              else
                relative_time_image_rgb(:, :, channel_index) = grayscale_image(:, :, 1);
              end
            end
            % whos relative_time_image_rgb

            relative_time_image_rgb_with_separate_images = relative_time_image_rgb;

            % Draw headings:

            draw_text_font_size = 1; pad_size = 15;
            draw_text_font_size = 3; pad_size = 20;
            top_heading_lines = 2;
            left_heading_lines = 3;

            % Allow for titles:
            top_heading_prefix = '';
            if ~isempty(video_titles)
                top_heading_lines = top_heading_lines + 1;
                top_heading_prefix = current_video_title;
                % relative_time_prefix = sprintf('Relative time %d', relative_time)
                % relative_time_prefix = sprintf('Relative time %+ 4d', all_relative_times_seconds_rounded(all_relative_time_index));
                relative_time_prefix = sprintf('Time %+ 4d', all_relative_times_seconds_rounded(all_relative_time_index));
                if ~isempty(top_heading_prefix)
                    top_heading_prefix = sprintf('%s; %s', relative_time_prefix, top_heading_prefix);
                else
                    top_heading_prefix = relative_time_prefix;
                end
            else
                % xruan remove cluster X in the title
                relative_time_prefix = sprintf('Time %+ 4d', timepoints_to_include(all_relative_time_index));
                top_heading_prefix = relative_time_prefix;
            end

            % Pad to make room for headings:
            relative_time_image_rgb_with_separate_images = padarray(relative_time_image_rgb_with_separate_images, pad_size * [top_heading_lines, left_heading_lines], 'pre');

            % Draw sensor headings:
            draw_text_options = struct('FontSize', draw_text_font_size);

            draw_text_centered_options = draw_text_options;
            draw_text_right_options = draw_text_options;
            draw_text_centered_options.alignment = 0;
            draw_text_right_options.alignment = 1;
            if isempty(video_titles)
                % Title:
                % beep, keyboard
                relative_time_image_rgb_with_separate_images(1:pad_size * 1, :, 1:3) = repmat(draw_text(top_heading_prefix, relative_time_image_rgb_with_separate_images(1:pad_size * 1, :, 1), [1 + pad_size * left_heading_lines + size(relative_time_image_rgb, 2) * 0, 1], draw_text_options), [1, 1, 3]);
                % Headings:
                relative_time_image_rgb_with_separate_images(:, :, 1) = draw_text({'', protein_names{1}}, relative_time_image_rgb_with_separate_images(:, :, 1), [1 + pad_size * left_heading_lines + size(relative_time_image_rgb, 2) * 0, 1], draw_text_options);
                % whos relative_time_image_rgb_with_separate_images
                if size(current_video_channels, 1) >= 2
                    relative_time_image_rgb_with_separate_images(:, :, 2) = draw_text({'', protein_names{2}}, relative_time_image_rgb_with_separate_images(:, :, 2), [1 + pad_size * left_heading_lines + floor(size(relative_time_image_rgb, 2) * .5), 1], draw_text_centered_options);
                end
                if size(current_video_channels, 1) >= 3
                    relative_time_image_rgb_with_separate_images(:, :, 3) = draw_text({'', protein_names{3}}, relative_time_image_rgb_with_separate_images(:, :, 3), [1 + pad_size * left_heading_lines + size(relative_time_image_rgb, 2) * 1, 1], draw_text_right_options);
                end
            else
                relative_time_image_rgb_with_separate_images(1:pad_size * 1, :, 1:3) = repmat(draw_text(top_heading_prefix, relative_time_image_rgb_with_separate_images(1:pad_size * 1, :, 1), [1 + pad_size * left_heading_lines + size(relative_time_image_rgb, 2) * 0, 1], draw_text_options), [1, 1, 3]);              
                relative_time_image_rgb_with_separate_images(:, :, 1) = draw_text(protein_names{1}, relative_time_image_rgb_with_separate_images(:, :, 1), [1 + pad_size * left_heading_lines + size(relative_time_image_rgb, 2) * 0, 1], draw_text_options);
                % whos relative_time_image_rgb_with_separate_images
                if size(current_video_channels, 1) >= 2
                    relative_time_image_rgb_with_separate_images(:, :, 2) = draw_text(protein_names{2}, relative_time_image_rgb_with_separate_images(:, :, 2), [1 + pad_size * left_heading_lines + floor(size(relative_time_image_rgb, 2) * .5), 1], draw_text_centered_options);
                end
                if size(current_video_channels, 1) >= 3
                    relative_time_image_rgb_with_separate_images(:, :, 3) = draw_text(protein_names{3}, relative_time_image_rgb_with_separate_images(:, :, 3), [1 + pad_size * left_heading_lines + size(relative_time_image_rgb, 2) * 1, 1], draw_text_right_options);
                end
            end

            if ~master_script_options.videos_color_headings
                relative_time_image_rgb_with_separate_images(1:pad_size * top_heading_lines, :, 2) = max(relative_time_image_rgb_with_separate_images(1:pad_size * top_heading_lines, :, :), [], 3);
                relative_time_image_rgb_with_separate_images(1:pad_size * top_heading_lines, :, 3) = max(relative_time_image_rgb_with_separate_images(1:pad_size * top_heading_lines, :, :), [], 3);
            end

            % Make the text more visible:
            relative_time_image_rgb_with_separate_images(1:pad_size * top_heading_lines, :, :) = relative_time_image_rgb_with_separate_images(1:pad_size * top_heading_lines, :, :) * 3;

            % Draw view headings (origin is below text, not above, and -3 is due to slight alignment issue):
            draw_text_options = struct('FontSize', draw_text_font_size, 'direction', 1);
            relative_time_image_rgb_with_separate_images(:, :, 1) = draw_text({'Slices', 'parallel to', 'top'}, relative_time_image_rgb_with_separate_images(:, :, 1), [1, 1 + pad_size * top_heading_lines + floor(size(relative_time_image_rgb, 1) * .5 * 1) - 3], draw_text_options);
            relative_time_image_rgb_with_separate_images(:, :, 1) = draw_text({'Slices', 'parallel to', 'synapse'}, relative_time_image_rgb_with_separate_images(:, :, 1), [1, 1 + pad_size * top_heading_lines + floor(size(relative_time_image_rgb, 1) * .5 * 2) - 3], draw_text_options);
            relative_time_image_rgb_with_separate_images(:, 1:pad_size * left_heading_lines, 2) = relative_time_image_rgb_with_separate_images(:, 1:pad_size * left_heading_lines, 1);
            relative_time_image_rgb_with_separate_images(:, 1:pad_size * left_heading_lines, 3) = relative_time_image_rgb_with_separate_images(:, 1:pad_size * left_heading_lines, 1);
            % Make the text more visible:
            relative_time_image_rgb_with_separate_images(:, 1:pad_size * left_heading_lines, :) = relative_time_image_rgb_with_separate_images(:, 1:pad_size * left_heading_lines, :) * 3;

            % Save a video frame for RGB and two-channel colorblind-friendly combinations in green and magenta as suggested at <http://jfly.iam.u-tokyo.ac.jp/color/index.html#stain>:

            % current_sensors = current_video_channels(:, 2);
            % % video_frame_filename = sprintf('%s_%s_%s_%s_frame%03d_%s', current_video_channels{1, 1}, current_sensors{:}, all_relative_time_index, model_type);
            % video_frame_filename = sprintf('%s_%s_%s_%s_frame%03d_%s_%s', current_video_channels{1, 1}, current_sensors{:}, all_relative_time_index, model_type, master_script_options.videos_contrast_enhancement_method);
            video_frame_filename = sprintf('%s_%s_frame%03d_%s', master_script_options.video_prefix, model_type, current_video_frame_index, master_script_options.videos_contrast_enhancement_method);
            video_frame_full_filename = [videos_location, video_frame_filename, '.png'];
            % imwrite(relative_time_image_rgb, video_frame_full_filename)
            imwrite(relative_time_image_rgb_with_separate_images, video_frame_full_filename)
        end
        
        video_filename = sprintf('%s_%s_%s.avi', master_script_options.video_prefix, model_type, master_script_options.videos_contrast_enhancement_method);

        % convert frames to a video 
         writerObj = VideoWriter([videos_location, video_filename]);
         writerObj.FrameRate = 1;
         % set the seconds per image
         % open the video writer
         open(writerObj);
        for current_video_frame_index = 1:num_timepoints_to_include
            video_frame_filename = sprintf('%s_%s_frame%03d_%s', master_script_options.video_prefix, model_type, current_video_frame_index, master_script_options.videos_contrast_enhancement_method);
            video_frame_full_filename = [videos_location, video_frame_filename, '.png'];
            cur_frame_image = imread(video_frame_full_filename);
            frame = im2frame(cur_frame_image);
            writeVideo(writerObj, frame);
         end
         % close the writer object
         close(writerObj); 
    end
end

end


function [model_relative_times, model_name_list] = extract_model_relative_times(model_name_list)

model_numbers = numel(model_name_list);
model_relative_times = zeros(model_numbers, 1);

for i = 1 : model_numbers
    cur_model_filename = model_name_list{i};
    try
        load(cur_model_filename);
    catch
        error('Unable to load model file!');
    end
    model_relative_times(i) = model.proteinModel.relative_time;
end
[model_relative_times, inds] = sort(model_relative_times);
model_name_list = model_name_list(inds);

end
 