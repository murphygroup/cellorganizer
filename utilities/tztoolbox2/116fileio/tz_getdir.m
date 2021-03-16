function dirpath = tz_getdir(dname)
%TZ_GETDIR
%   DIRPATH = TZ_GETDIR(DNAME)
%   
%   See also

%   11-Jul-2006 Initial write T. Zhao
%   Copyright (c) 2006 Murphy Lab
%   Carnegie Mellon University
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation; either version 2 of the License,
%   or (at your option) any later version.
%   
%   This program is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   General Public License for more details.
%   
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
%   02110-1301, USA.
%   
%   For additional information visit http://murphylab.web.cmu.edu or
%   send email to murphy@cmu.edu


if nargin < 1
    error('Exactly 1 argument is required');
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
dirpath = [];

[dname,dname2] = strtok(dname,'/');

switch dname
    case {'run','data'}
        dirpath = [sharedPath '/' dname];
    case 'fun'
        dirpath = [sharedPath '/' glFundir];
    case 'const'
        dirpath = [sharedPath '/' glConstdir];
    case 'cvs'
        dirpath = slicpath;
    case 'lib'
        dirpath = [sharedPath '/' glDowndir];
    case 'docs'
        dirpath = [sharedPath '/' glDocdir];
    case 'mt'
        dirpath = glMtpath;;
    case 'shared'
        dirpath = sharedPath;
    case 'local'
        dirpath = [glMtpath filesep glMachine];
    case 'testing'
        dirpath = [sharedPath '/run/' dname];
    otherwise
        switch dname(1)
            case '*'
                dirpath = [slicpath '/' dname(2:end) '/' 'matlab'];
            case '~'
                dirpath = [glMtpath filesep glMachine filesep dname(2:end)];
            otherwise
                dirname=[sharedPath '/' glFundir '/116' dname];
                if ~exist(dirname,'dir')                   
                    functionPath = which(dname);
                    dirname = fileparts(functionPath);
                    
                    if ~exist(dirname,'dir')   
                        dirpath = [];
                        return;
                    end
                end

                dirpath = dirname;
        end
end

if ~isempty(dname2)
    dirpath = [dirpath dname2];
end
