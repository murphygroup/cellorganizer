function list = tz_ls(pattern,option)
%TZ_LS List files or directories.
%   LIST = TZ_LS(PATTERN) list all files and directories matching PATTERN.
%   The results are stored in the cell array LIST.
%   
%   LIST = TZ_LS(PATTERN,OPTION) lets users specify listing files or 
%   directories by OPTION:
%       'all' - list all
%       'file' - list files
%       'dir' - list directories
%       '-r' - go through subdirectories
%   These option can be combined into one string, but each option unit
%   should be separated by a space. For example, a valid option could 
%   be 'all -r', which means list all files and directories
%   without dot under the directory and its subdirectories provided in
%   pattern.
%   A string in LIST is a fullpath only when option has '-r' and the
%   corresponding file or directory is not under the current directory
%   directly.
%
%   NOTICE: current version only support following combinations:
%   'all -r', 'file -r', 'dir -r'

%   11-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('1 or 2 arguments are required')
end

if nargin < 2
    option = 'all';
end

options = tz_strtok(option,' ');

subDirOptionIndex = strmatch('-r',options,'exact');
if(isempty(subDirOptionIndex))
    patterns{1} = pattern;
else
    options{subDirOptionIndex}=[];
        
    if strmatch('./',pattern)
        currentDir = pwd;
        if currentDir(end) ~= filesep
            currentDir(end+1) = filesep;
        end
        pattern = strrep(pattern,'./',currentDir);
    end
    
    [directory,file,extension] = fileparts(pattern);
    filePattern = [file,extension];
    
    if exist(fullfile(directory,filePattern),'dir')
        directory = fullfile(directory,filePattern);
        filePattern = '';
    end
    
    isAbsolutePath = 1;
    if isempty(directory)
        isAbsolutePath = 0;
    else
        if directory(1) ~= filesep
            isAbsolutePath = 0;
        end
    end
    
    if ~isAbsolutePath
        currentDir = pwd;
        directory = fullfile(currentDir,directory);
    end
    
    dirs={};
    dirs = tz_lsdir(directory,dirs);
    dirs = {directory,dirs{1:end}};
    
    for i=1:length(dirs)
        patterns{i} = fullfile(dirs{i},filePattern);
    end
end 
  
list = {};
  
if isempty(options{1})
    options{1} = 'all';
end

for iPattern=1:length(patterns)
    s = dir(patterns{iPattern});
    
    for i=1:length(s)       
        if exist(patterns{iPattern},'dir')
            currentDir = patterns{iPattern};
        else
            currentDir = fileparts(patterns{iPattern});
        end
        
        if strcmp(currentDir,pwd)
            currentDir = '';
        end
        
        if strcmp(s(i).name,'.') | strcmp(s(i).name,'..') 
            continue
        end
        
        if(isempty(subDirOptionIndex))
            switch options{1}
            case 'all'
                list{end+1} = s(i).name;
            case 'file'
                if ~s(i).isdir
                    list{end+1} = s(i).name;
                end
            case 'dir'
                if s(i).isdir
                    list{end+1} = s(i).name;
                end
            end       
            
        else
            switch options{1}
            case 'all'
                list{end+1} = fullfile(currentDir,s(i).name);
            case 'file'
                if ~s(i).isdir
                    list{end+1} = fullfile(currentDir,s(i).name);
                end
            case 'dir'
                if s(i).isdir
                    list{end+1} = fullfile(currentDir,s(i).name);
                end
            end
        end
%         case 'nodot'
%             if s(i).name(1) ~= '.'
%                 list{end+1} = s(i).name;
%             end 
%         end
    end
end

list = list';


%%%%%%%%%%%%%%sub function%%%%%%%%%%%%%%%%%%
function list = tz_lsdir(directory,list)
%TZ_LSDIR List all sub directories.

ss = dir(directory);
hasMessed = 1;

for i=1:length(ss)
    if hasMessed
        if strcmp(ss(i).name,'.')
            ss(i).isdir = 0;
        end
        if strcmp(ss(i).name,'..') 
            ss(i).isdir = 0;
            hasMessed = 0;
        end
    end
    
    if ss(i).isdir
        list{end+1} = fullfile(directory,ss(i).name);
        list = tz_lsdir(list{end},list);
    end
end
