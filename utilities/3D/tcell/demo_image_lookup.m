function demo_image_lookup
% the demo is used for check the video frames. 
% one can choose specify frame and specify slices to generate a video 

% Author: Xiongtao Ruan xruan@andrew.cmu.edu
% Date: Nov. 16, 2015


[curr_dir] = fileparts(mfilename('fullpath'));
cd(curr_dir);

dir_info = dir_recursive(curr_dir);

dir_names = unique(cellfun(@(x) fileparts(x), dir_info, 'uniformoutput', false))';

dir_names(cellfun(@(x) strcmp(x, curr_dir), dir_names)) = [];
dir_names = dir_names(~cellfun(@isempty, regexpi(dir_names, 'GFP')));

disp(sprintf('There are %d sets of video(s), the names are:', numel(dir_names)))
for i = 1 : numel(dir_names);
    disp(sprintf('%d: %s', i, dir_names{i}));
end

num = input('You can choose an interested video, and please type in a number which stands for the video above\n');

disp(sprintf('We are going to play with video %d', num));

disp_dir = dir_names{num};

disp_frames = dir_recursive(disp_dir);

N_frame = numel(disp_frames);

frame_timepoints = regexpi(disp_frames, 'Timepoint ?([0-9]+)|T([0-9]+)', 'tokens');
frame_timepoints = frame_timepoints(~cellfun(@isempty, frame_timepoints));
frame_timepoints = cell2mat(cellfun(@(x) str2num(x{1}{1}), frame_timepoints, 'uniformoutput', false));
disp(sprintf('There are totally %d frames', numel(frame_timepoints)));

% disp('Do you want a specific frame or choose a specific slice from some frames?');
% 
% mode = input('please type F: frame, S: slice', 's');
% 
% if strcmpi(mode, 'f')

% get frame numbers
frame_num = input('please specify which frame(s) do you want? use "," to separate frame number, \n and use "-" to represent the sequence, eg 1-9 stands for 1, 2, ..., 9\n', 's')
frame_num = strsplit(frame_num, ',');
frame_num_int = [];
for i = 1 : numel(frame_num);
    num_i = frame_num{i};
    if isempty(regexp(num_i, '-'))
        frame_num_int = [frame_num_int, str2num(num_i)];
    else
        num_se = strsplit(num_i, '-');
        frame_num_int = [frame_num_int, str2num(num_se{1}) :  str2num(num_se{2})];
    end       
end

frame_num_int = sort(unique(frame_num_int));

% get slice numbers
slice_num = input('please specify which slice(s) do you want? use "," to separate slice number, \n and use "-" to represent the sequence, eg 1-9 stands for 1, 2, ..., 9, and use "ALL" to display all slices\n', 's')

slice_num_int = [];

if regexpi(slice_num, 'all')
    slice_num_int = -1;
else
    slice_num = strsplit(slice_num, ',');
    slice_num_int = [];
    for i = 1 : numel(slice_num);
        num_i = slice_num{i};
        if isempty(regexp(num_i, '-'))
            slice_num_int = [slice_num_int, str2num(num_i)];
        else
            num_se = strsplit(num_i, '-');
            slice_num_int = [slice_num_int, str2num(num_se{1}) :  str2num(num_se{2})];
        end       
    end
end

slice_num_int = sort(unique(slice_num_int));


image_sets = dir_recursive(disp_dir);
image_used = {};
for i = 1 : numel(frame_num_int)
    current_frame = frame_num_int(i);
    if min(frame_timepoints) == 0
        current_frame = current_frame - 1;
    end
    current_frame_name = image_sets(~cellfun(@isempty, regexpi(image_sets, sprintf('Timepoint ?%02d|T%02d', current_frame, current_frame))));
    image_i = current_frame_name{1};
    
    im_info = imfinfo(image_i);
    num_img = numel(im_info);

    for j = 1 : num_img
        I_j = double(imread(image_i, j));
        I(:, :, j) = double(I_j);
    end
    image_used{i} = mat2gray(I);
end

% save the video 
now_time = datestr(now, 'yyyyMMDD_HHmmss');
video_filename = sprintf('%s/video_%s.avi', curr_dir, now_time);

if slice_num_int ~= -1
    image_used = cellfun(@(x) x(:, :, slice_num_int), image_used, 'uniformoutput', false);
end
create_video(video_filename, image_used)    

end


function create_video(filename, image_sets)

v = VideoWriter(filename,'Uncompressed AVI');
open(v);
for i = 1 : numel(image_sets)
    for j = 1 : size(image_sets{i}, 3)
        writeVideo(v, image_sets{i}(:, :, j))
    end
end
close(v);
end

function [filenames] = dir_recursive(given_dir)

    % xruan 08/13/2015
    % get all filenames under the current dir and its subdiretories 
    % use DFS to recursively get all filenames (no directory)
    
    if strcmp(given_dir(end), filesep)
        given_dir(end) = [];
    end
      
    filenames = {};
    
    stack = {};
    
    if isdir(given_dir)
        stack{end + 1} = given_dir;
    end
    
    
    while ~isempty(stack)

        % pop one element
        current_dir = stack{end};
        stack(end) = [];

        fileinfo = dir(current_dir);
        if isempty(fileinfo)
            continue;
        end
        subdir_files = {fileinfo.name};
        sub_isdir = [fileinfo.isdir];
        subdir_files(1 : 2) = [];
        sub_isdir(1 : 2) = [];
        subdirs = cellfun(@(x) [current_dir, filesep, x], subdir_files(sub_isdir), 'UniformOutput', false);
        if length(subdirs) > 0
            stack(end + 1 : end + length(subdirs)) = subdirs;
        end

        subfiles = cellfun(@(x) [current_dir, filesep, x], subdir_files(~sub_isdir), 'UniformOutput', false);
        if length(subfiles) > 0
            filenames(end + 1 : end + length(subfiles)) = subfiles;
        end
    end
end



















