function gd(dname)
%GD Shortcut of change current directory.
%   GD(DNAME) has the following shortcuts with DNAME:
%       '-b': back
%       '-f': forward
%       '-l': list directories in history
%       'mt': /home/tingz/matlab
%       'cvs': /home/tingz/matlab/newcvs/NEWSLIC
%       'fun': /home/tingz/matlab/toolbox
%       'run': /home/tingz/matlab/run
%       'const': /home/tingz/matlab/constant
%       'docs': /home/tingz/matlab/doc
%       'lib': /home/tingz/matlab/download
%   If DNAME starts with '*', gd will go to cvs subdirectory. If it starts
%   with '~', gd will go to local subdirectory. For example:
%       gd *classify - go to /home/tingz/matlab/cvs/NEWSLIC/classify/matlab
%       gd ~run - go to /home/tingz/matlab/goblin/run (on goblin)
%   
%   If DNAME does not match all directories described above, it will be
%   take as a [lib] directory. But if there is still no such a directory,
%   DNAME will then be interpreted as a function name and look for in in
%   current paths. If the function is found, it will go to the directory 
%   where the function is located. For example:
%       gd classif - go to /home/tingz/matlab/shared/toolbox/116classif
%       gd tz_exp - go to /home/tingz/matlab/shared/toolbox/116math
%
%   See also HD

%   06-Apr-2005 Initial write T. Zhao
%   03-Mar-2006 Modified T. Zhao
%       - support function name
%       - remove appending space at the end
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argumentis required')
end

global glMtpath glShareddir glMachine glCvsdir glFundir glDatadir ...
    glRundir glConstdir glDocdir glDowndir glDirhistprev glDirhistnext;

% tz- 03-Mar-2006
% if length(dname)<3
%     dname=[dname ' '];
% end
% tz--

sharedPath = [glMtpath filesep glShareddir];

slicpath=[glMtpath filesep glCvsdir filesep 'SLIC'];

switch dname
    case '-b'
        if ~isempty(glDirhistprev)
            glDirhistnext{end+1} = pwd;
            cd(glDirhistprev{end});
            glDirhistprev(end) = [];
        else
            disp('No back directory!');
        end
    case '-f'
        if ~isempty(glDirhistnext)
            glDirhistprev{end+1} = pwd;
            cd(glDirhistnext{end});
            glDirhistnext(end) = [];
        else
            disp('No forward directory!');
        end
    case '-l'
        for i=1:length(glDirhistprev)
            disp(glDirhistprev{i});
        end
        disp(['--> ' pwd]);
        for i=1:length(glDirhistnext)
            disp(glDirhistnext{i});
        end
%         disp(glDirhistprev')
%         disp(glDirhistnext')
    otherwise
        hd(tz_getdir(dname));
end

