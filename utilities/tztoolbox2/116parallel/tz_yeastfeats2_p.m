function tz_yeastfeats2_p(datafile)
%TZ_YEASTFEATS_P Calculate features for yeast images.
%   TZ_YEASTFEATS_P(DATAFILE) saves the calcuate features under the
%   directory specified in the data file DATAFILE, which is a mat file
%   containing two variables:
%       'resultdir' - directory for saving calculated features
%       'procfiles' - image files to calcualte features
%   
%   See also

%   22-Nov-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('Exactly 1 arguments are required')
end

if license('checkout','statistics_toolbox')==0
  return
end

addpath /home/tingz/matlab/shared
tz_initpath

load(datafile)
nprocfiles = length(procfiles{1})
isbreak = 0;

for n=1:nprocfiles
    imagePath = strrep(procfiles{3}{n},'matlab/data','tmp');
    featureFileName = tz_yeastimg2featname(imagePath);
    controlDirectory = [resultdir filesep featureFileName];
    [s,msg] = mkdir(resultdir,featureFileName);
        
    %If the job has been taken over
    if strfind(msg,'exist')
        continue
    end
       
    n
    %%%%%%%%TO DO%%%%%%%%%%
    resultfile = [resultdir filesep ...
        tz_yeastimg2featname(imagePath) '.mat'];
    if ~exist(resultfile,'file');
%         features = n;
        img = ml_readimage(imagePath);
        dnaImagePath = strrep(procfiles{1}{n},'matlab/data','tmp');
        dnaimg = ml_readimage(dnaImagePath);
        [names, features, slfnames] = ml_featset( double(img), [], ...
            double(dnaimg), 'all180',0,0,'quantile',...
            'rc');
        
        save(resultfile,'features');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%
    rmdir(controlDirectory)
end
