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


