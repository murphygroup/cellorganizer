function tz_init_yeastfeats_p(yeastRoot)
%TZ_INIT_YEASTFEATS_P Initialize data for yeast image analysis.
%   TZ_INIT_YEASTFEATS_P(YEASTROOT) creates several files for yeast feature
%   calculation.
%   
%   See also TZ_YEASTFEATS_P

%   23-Nov-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

% yeastRoot = '/home/tingz/matlab/data/yeastGFP';
addpath /home/tingz/matlab
tz_initpath

cloneDirs = tz_ls(yeastRoot,'dir');
nimage = [];
k=1;
procfiles{1} = {};
procfiles{2} = {};
procfiles{3} = {};
% for i=1:length(cloneDirs)
% imgdir = [yeastRoot filesep cloneDirs{i}];
procfiles{1} = tz_ls([yeastRoot filesep '*_DAPI.png'],'-r');
procfiles{2} = tz_ls([yeastRoot filesep '*_DIC.png'],'-r');
procfiles{3} = tz_ls([yeastRoot filesep '*_GFP.png'],'-r');

% procfiles{1} = {procfiles{1}{:} imgFiles{i}{1}{:}};
% procfiles{2} = {procfiles{2}{:} imgFiles{i}{2}{:}};
% procfiles{3} = {procfiles{3}{:} imgFiles{i}{3}{:}};
% end

resultdir = '/home/tingz/matlab/data/yeastfeats';

ctrfile = [resultdir filesep 'tz_yeastfeats_p.ctr'];
prvfile = [ctrfile '.prv'];

fid = fopen(ctrfile,'w');
fwrite(fid,uint16(0),'uint16');
fclose(fid)

datafile = '/home/tingz/matlab/data/procfiles.mat';
save(datafile, ...
    'procfiles','resultdir','ctrfile','prvfile','yeastRoot');

