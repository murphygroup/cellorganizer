function [status] = mkdir_recursive(full_path)
    % xruan 07/25/2015 
    % mkdir recursively
    path_list = strsplit(full_path, '/');

    path_depth_i = '';
    for i = 1 : length(path_list)
        if strcmp(path_list{i}, '~')
            path_depth_i = '~';
            continue;
        end
        path_depth_i = strjoin( {path_depth_i, path_list{i}}, '/');
        if ~exist(path_depth_i, 'dir')
            mkdir(path_depth_i);
        end

    end
    status = 1;

end