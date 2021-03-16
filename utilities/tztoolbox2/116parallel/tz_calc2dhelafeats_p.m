function tz_calc2dhelafeats_p(paramfile)
%TZ_CALC2DHELAFEATS_P Calculate features for 2D Hela images
%   TZ_CALC2DHELAFEATS_P(PARAMFILE)
%   
%   See also

%   18-May-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('Exactly 1 argument is required')
end

if license('checkout','statistics_toolbox')==0
  return
end

addpath /home/tingz/matlab/shared
tz_initpath

param = load(paramfile);

classes = tz_ls(param.imgdir,'dir');
nClass = length(classes);

channels = {'crop','dna','prot'};
resultdir = param.resultdir;

for iClass = 1:nClass

    classDir= [param.imgdir filesep classes{iClass}];
    proteinDir = [classDir filesep channels{3}];
    cropDir = [classDir filesep channels{1}];
    dnaDir  =  [classDir filesep channels{2}];
      
    proteinFiles = tz_ls([proteinDir filesep '*.dat']);
    cropFiles = tz_ls([cropDir filesep '*.tif*']);
    dnaFiles = tz_ls([dnaDir filesep '*.dat']);
    nCell = length(proteinFiles);
    features = [];
    for iCell = 1:nCell   
        featfile = [classes{iClass} num2str(iCell) '.mat'];
        controlDirectory = [resultdir filesep featfile '.ctr'];
        [s,msg] = mkdir(resultdir,[featfile '.ctr']);
        %If the job has been taken over
        if strfind(msg,'exist')
            continue
        end
        resultFile = [resultdir filesep featfile]
        
        if ~exist(resultFile,'file')
            proteinPath = [proteinDir filesep proteinFiles{iCell}];
            cropPath = [cropDir filesep cropFiles{iCell}];
            dnaPath = [dnaDir filesep dnaFiles{iCell}];
            proteinImage = ml_readimage(proteinPath);
            cropImage = ml_readimage(cropPath);
            dnaImage =  ml_readimage(dnaPath);
            [feat_names,features,feat_slf] = ...
                ml_featset(proteinImage, ...
                cropImage, dnaImage,param.featset);
            tz_save(resultFile,{feat_names,features,feat_slf}, ...
                {'feat_names','features','feat_slf'}, ...
                ['tz_calc2dhelafeats_p(''' paramfile ''')'], ...
                ['features for ' proteinPath]);
        end
        rmdir(controlDirectory);
    end
    
    allfeatures{iClass} = features;
end
