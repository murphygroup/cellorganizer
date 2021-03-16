function tz_calcgenimgfeats(paramfile)
%TZ_CALCGENIMGFEATS Calculation features for generated images
%   TZ_CALCGENIMGFEATS(PARAMFILE) save feature files in a directory
%   specified in PARAMFILE, which is a MAT file.
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

imageFiles = tz_ls([param.imgdir filesep param.pattern]);
resultdir = param.resultdir;

for i=1:length(imageFiles)
    imagePath = [param.imgdir filesep imageFiles{i}];
    featfile = [imageFiles{i} '.mat'];
    controlDirectory = [resultdir filesep featfile '.ctr'];
    [s,msg] = mkdir(resultdir,[featfile '.ctr']);
    %If the job has been taken over
    if strfind(msg,'exist')
        continue
    end
    resultFile = [resultdir filesep featfile]
    
    if ~exist(resultFile,'file')
        %calculate features
        img = ml_readimage(imagePath);
        
        [feat_names, features, feat_slf] = ...
            ml_features(img*255, [], img>0,...
            {'skl' 'img' 'hul' 'zer' 'har' 'edg'});
        tz_save(resultFile,{feat_names,features,feat_slf}, ...
            {'feat_names','features','feat_slf'}, ...
            ['tz_calcgenimgfeats(''' paramfile ''')'], ...
            ['features for ' imagePath]);
    end
    rmdir(controlDirectory);
end

