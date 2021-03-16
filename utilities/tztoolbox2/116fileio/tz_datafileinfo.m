function info = tz_datafileinfo(filepath,option)
%TZ_DATAFILEINFO Show the information of a data file.
%   INFO = TZ_DATAFILEINFO(FILEPATH) returns the information of the data file
%   FILEPATH and displays the information on the screen as well.
%   
%   INFO = TZ_DATAFILEINFO(FILEPATH,OPTION) selects which informatin to show
%   according to the string OPTION:
%       'script' - the script that created the file
%       'comments' - the author's comments about the file
%       'machine' - the machine in which the file is created
%       'savetime' - the time when the file is created
%       'matlab' - the version of which the Matlab is used to create the file
%       'version' - the version of the file
%       'summary' - a summary of the data file
%
%   See also

%   26-Oct-2006 Initial write T. Zhao
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

data = load(filepath);

switch option
    case {'script','comments','machine','savetime','matlab','version'}
        if isfield(data,option)
            info = getfield(data,option);
            if strcmp(option,'savetime')
                disp(datestr(info));
            else
                disp(info);
            end
            
        else
            info = [];
            disp([option ' is not available']);
        end
    case 'summary'
        disp(data);
    otherwise
        error(['Unrecognized option: ' option]);
end

