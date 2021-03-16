function tz_init_yeastfeats_p(yeastRoot,resultdir)
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


% resultdir = '/home/tingz/matlab/data/yeastfeats';
unix(['mkdir ' resultdir]);
unix(['mkdir ' resultdir filesep 'seg']);

datafile = [resultdir filesep 'procfiles.mat'];
save(datafile, ...
    'procfiles','resultdir','yeastRoot');

