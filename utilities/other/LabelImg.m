function LabelImg(filePath,varname,filetype)
%Searches all folders and subfolders, reads mat files and plots each
%variable in the mat file and asks the user for a label for the data
%Inputs: 
%filePath = string pointing to top directory you want to start searching
%through
%varname = name of variable from mat to read/plot, if none specified will
%read and plot all files in whatever mat file are found 
%filetype = type of file to read. can be 'mat' or 'tif' currently


if nargin <3
    filetype = 'mat';
end

%first get the files that end in .mat from the current path
[s,filesstr] = system(['find ',filePath,' -type f -name "*.' filetype '"']);

%loop through all the files and separate them from one char array into
%a cell array of strings for reading
files = textscan(filesstr,'%s');
files = files{1};

filelabels = 1;
filenames = cell(1,length(files))
for i = 1:length(files)
    if isempty(varname)
        if strcmp(filetype,'mat')
            currimgs = load(files{i});
            imgnames = fieldnames(currimg);
            numoffiles = length(imgnames);
        elseif strcmp(filetype,'tif')
            currimgs = tif2img(files{i});
            imgnames = [];
            numoffiles = 1;
        end
    else
        currimgs = load(files{i},varname);
        imgnames = fieldnames(currimg);
        numoffiles = length(imgnames);
    end
    
    for j = 1:length(numoffiles)
        if isstruct(currimgs)
            currimg = getfield(currimgs,imgnames{j});
        else
            currimg = currimgs;
        end
%         h = plot(currimg);
        h = im2projection(currimg);
        label{k} = input('Is it a good segmentation? Y/N:','s');
        currfile{k} = [files,imgnames{j}];
        k = k+1;
        close(h)
        
    end
    
end
