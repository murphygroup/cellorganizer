function [out_stat, deleted_files] = chunk_lock_clean(path, duration)
% This function is for cleaning tmp files whose jobs are not running in the folder.
%

% Author: Xiongtao Ruan
% Date: Feb. 24, 2016

if nargin < 2
    delete_by_time = false;
else
    delete_by_time = true;
end

deleted_files = {};  
out_stat = false;

[status, output] = system(sprintf('find %s -type f -name "*.tmp"', path));
if ~isempty(output)
    temp_filenames = strsplit(strtrim(output), '\n');
    if numel(temp_filenames) > 0 
        for i = 1 : numel(temp_filenames)
            current_temp_filename = temp_filenames{i};
            if delete_by_time
                temp_file_info = dir(current_temp_filename);
                if (datenum(datetime('now')) - [temp_file_info.datenum]) * 24 * 60 > duration
                    deleted_files{end + 1} = current_temp_filename;
                    delete(current_temp_filename);
                end
            else
                % current_temp_filename
                try
                    job_info = textread(current_temp_filename, '%s', 'delimiter', '\n');
                catch
                    if ~exist(current_temp_filename, 'file')
                        continue;
                    end
                end
                if numel(job_info) ~=2
                    deleted_files{end + 1} = current_temp_filename;
                    delete(current_temp_filename);
                    continue
                end
                id_type = job_info{1};
                id = job_info{2};
                 
                % for the cluster to check job id
                if strcmp(id_type, 'jid')
                    [status, output] = system(sprintf('sacct -j %s --format=State', id));           
                
                    if isempty(regexp(output, 'RUNNING'))  
                        deleted_files{end + 1} = current_temp_filename;
                        delete(current_temp_filename);
                    end
                elseif strcmp(id_type, 'pid')
                    [status, output] = system(sprintf('ps -o pid= -p %s', id));  
                    if isempty(output)
                        deleted_files{end + 1} = current_temp_filename;
                        delete(current_temp_filename);
                    end
                else
                    deleted_files{end + 1} = current_temp_filename;
                    delete(current_temp_filename);
                end
            end
        end
    end
end

if ~isempty(deleted_files)
    out_stat = true;
end

end
    

