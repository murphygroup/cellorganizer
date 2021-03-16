function [stat,msg] = ml_mkdir(dirname,hr)
%ML_MKDIR Make directory.
%   STAT = ML_MKDIR(DIRNAME) tries to make directory DIRNAME and returns the
%   status of the result.
%   
%   STAT = ML_MKDIR(DIRNAME,HR) makes directories hierachically if HR is true.
%   Otherwise it is the same as ML_MKDIR(DIRNAME).
%   
%   [STAT,MSG] = ML_MKDIR(...) also returns the message after the attepmt.
%   
%   See also

%   30-Jan-2007 Initial write T. Zhao
%   Copyright (c) 2007 Murphy Lab
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

if ~exist('hr','var')
    hr = 0;
end

if ~hr
    [stat,msg] = mkdir(dirname);
else
    dirs = {};
    if dirname(end) == filesep
        dirname(end) = [];
    end
    if ~isempty(dirname)
        while ~exist(dirname,'dir') 
            [dirname curdir] = fileparts(dirname);
            dirs = {dirs{:} curdir};
        end
        if ~isempty(dirs)
            for i=1:length(dirs)
                dirname = fullfile(dirname,dirs{i});
                [stat,msg] = mkdir(dirname);
            end
        else
            stat = 0;
            msg = 'The directory alread exists';
        end
    else
        stat = 0;
        msg = 'invalid directory name';
    end
end
