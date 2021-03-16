function tz_help
%TZ_HELP Open help file of T. Zhao's toolbox.
%   TZ_HELP is provides a convenient way to learn how to use T. Zhao's
%   toolbox.

%   27-Aug-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

global glMtpath glDocdir glShareddir

helpFile = [glMtpath filesep glShareddir filesep glDocdir '/index.html'];

disp('Help starting ...');
[s,msg] = unix('which mozilla');

if s==0
    tz_open(helpFile);
else
    open(helpFile);
end
