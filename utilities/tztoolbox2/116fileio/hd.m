function hd(pathname)
%HD Change directory and record.
%   HD(PATHNAME) does the same thing as CD. But the directory changing
%   history will be recorded.
%   
%   See also GD

%   15-Aug-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argument is required')
end

global glDirhistprev glDirhistnext;

if exist(pathname,'dir')
    glDirhistprev{end+1} = pwd;
    glDirhistnext = {};
    cd(pathname);
else
    disp('no such directory');
end