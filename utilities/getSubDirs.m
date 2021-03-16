function [ folderPaths, folderNames] = getSubDirs( parentDir, wildcard )
%Gregory Johnson - 9/20/2011

warning off;

searchDir = parentDir;

if exist('wildcard', 'var')
    searchDir = [searchDir filesep wildcard];
end
    

topFolders = dir(searchDir);
counter = 1;

%remove all non-dir files

topFolders = topFolders([topFolders(:).isdir]);

if ~exist('wildcard', 'var')
    topFolders = topFolders(3:end);
end

folderNames = {topFolders(:).name};

folderPaths = cell(length(folderNames),1);

for i = 1:length(folderNames)
    folderPaths{i} = [parentDir filesep folderNames{i}];
end


end

