function tz_maketoolboxdoc
%TZ_MAKETOOLBOXDOC makes documents for the toolbox
%   TZ_MAKETOOLBOXDOC processes files in T. Zhao's toolbox and creates help
%   files.
%   
%   See also m2html

%   07-Aug-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

global glMtpath glShareddir glFundir glDocdir

prevpath=pwd;
glSharedpath = [glMtpath filesep glShareddir];

cd(glSharedpath)
rmpath(tz_getdir('~toolbox'));

m2html('mfiles',glFundir,'htmldir',[glDocdir '/' glFundir], ...
    'recursive','on','source','off');
tz_prochelpindex([glSharedpath '/' glDocdir '/' glFundir ...
    '/' 'index.html'])

addpath(tz_getdir('~toolbox'));

cd(prevpath)
