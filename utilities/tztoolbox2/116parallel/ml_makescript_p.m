function ml_makescript_p(script,dirname)
%ML_MAKESCRIPT_P Create a script for parallel processing.
%   ML_MAKESCRIPT_P(SCRIPT) create a script for parallel processing. SCRIPT is 
%   a script name for running. It must be a script that can be run anywhere
%   after the toolbox is initialized. The new script will be created under
%   ~toolbox.
%   
%   ML_MAKESCRIPT_P(SCRIPT,DIRNAME) creates a script under the directory 
%   DIRNAME.
%
%   The grammer of making a script for parallel processing:
%       '%@parallel'
%       '@control'
%       '@release'
%   
%   See also

%   21-Aug-2006 Initial write T. Zhao
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
    error('1 or 2 arguments are required');
end

scriptPath = which(script);

if ~exist(scriptPath,'file')
    disp(['The script ''' script ''' can not be found']);
    return;
end

if ~exist('dirname','var')
    dirname = tz_getdir('~toolbox');
end

pscriptPath = [dirname filesep script '_p.m'];
fid = fopen(pscriptPath,'w');
codeLine = ['addpath ' tz_getdir('shared')]; 
fprintf(fid,'%s\n',codeLine);
codeLine = 'tz_initpath';
fprintf(fid,'%s\n',codeLine);

fid2=fopen(scriptPath);
while 1
    tline = fgetl(fid2);
    controlIndex = strfind(tline,'%@parallel'); %label of parallel command
    if ~isempty(controlIndex)
        tokens = tz_strtok(tline(controlIndex:end),' ');
        indent = tline(1:controlIndex-1);
        switch tokens{2}
            case '@control' %label of secondary command for control
                controlDirectory = ['[' tokens{3} ' ''.ctr'']'];
                fprintf(fid,'%s\n', ...
                    [indent '[s,msg] = mkdir(' controlDirectory ');']);
                fprintf(fid,'%s\n',[indent 'if strfind(msg,''exist'')']);
                fprintf(fid,'%s\n',[indent '    continue;']);
                fprintf(fid,'%s\n',[indent 'end']);
            case '@release' %label of secondary command for release
                if exist('controlDirectory','var')
                    fprintf(fid,'%s\n', ...
                        [indent 'rmdir(' controlDirectory ')']);
                end
            otherwise %additional command lines
                fprintf(fid,'%s\n', ...
                    [indent tz_cell2str(tokens(2:end),' ')]);
        end 
    else %original codes
        if ~ischar(tline), break, end
        fprintf(fid,'%s\n',tline);
    end
end
fclose(fid2);

% fprintf(fid,'%s\n',script);

fclose(fid);

disp(['The script ' pscriptPath ' is created successfully from ' ...
      scriptPath]);

