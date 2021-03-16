function [ files ] = loadFiles( file_paths, variables )
%LOADFILES Summary of this function goes here
%   Detailed explanation goes here

if ~exist('variables', 'var')
    variables = '';
end

if ischar(variables)
    variables = {variables};
end
    
files = cell(1,length(file_paths));

for i = 1:length(file_paths)
    if exist(file_paths{i}, 'file')
        out_tmp = load(file_paths{i}, variables{:});
        
        if ~isempty(out_tmp)
            files{i} = out_tmp;
        end
    end
end

end

