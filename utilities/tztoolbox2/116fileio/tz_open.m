function tz_open(filename)
%TZ_OPEN Open a file.
%   TZ_OPEN(FILENAME) opens a file with the name FIELNAME. If the file has
%   extension '.htm' or '.html', it will be opened in Mozilla. Otherwise,
%   it is opened by emacs.

%   13-DEC-2004 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argumentis required')
end

fullpath=which(filename);
if isempty(fullpath)
    disp(['Script ' filename ' not found']);
else
    [pathstr,name,ext] = fileparts(fullpath);
    if strcmp('.htm',ext) | strcmp('.html',ext)
        openSoftware = 'mozilla';
    else
        openSoftware = 'emacs';
    end
    cmd = [openSoftware ' ' fullpath ' &'];
    unix(cmd);
    disp(['Opening ' fullpath ' ...']);
end
