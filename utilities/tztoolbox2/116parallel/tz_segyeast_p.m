function tz_segyeast_p(datafile)
%TZ_YEASTFEATS_P A script for segmenting yeast images in clusters.
%   TZ_YEASTFEATS_P(DATAFILE,RESULTDIR)
%   
%   See also

%   22-Nov-2005 Initial write T. Zhao
%   16-Feb-2005 Modified T. Zhao
%       - change '_seg' into '.seg'
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('Exactly 1 arguments are required')
end

if license('checkout','statistics_toolbox')==0
  return
end

addpath /home/tingz/matlab
tz_initpath

load(datafile)
nprocfiles = length(procfiles{1})
isbreak = 0;

for n=1:nprocfiles
    dicImagePath = strrep(procfiles{2}{n},'matlab/data','tmp');
    segFileName = [tz_yeastimg2featname(dicImagePath) '.seg'];
    controlDirectory = [resultdir filesep 'seg' filesep segFileName];
    [s,msg] = mkdir([resultdir filesep 'seg'],segFileName);
        
    %If the job has been taken over
    if strfind(msg,'exist')
        continue
    end
       
    n
    %%%%%%%%TO DO%%%%%%%%%%
    resultfile = [resultdir filesep 'seg' filesep...
        segFileName '.mat'];
    if ~exist(resultfile,'file')
        dnaImagePath = strrep(procfiles{1}{n},'matlab/data','tmp');
        segimg = uint16(ml_cutcells(dnaImagePath,dicImagePath));

        save(resultfile,'segimg');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%
    unix(['rm -r ' controlDirectory]);
end
