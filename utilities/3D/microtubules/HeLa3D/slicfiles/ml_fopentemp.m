function [arg1, arg2] = ml_fopentemp( dirname)

% [FILEID, FILENAME] = ML_FOPENTEMP( DIRNAME)
% [FILENAME] = ML_FOPENTEMP( DIRNAME)
%
% Creates and Opens a temporary file in directory DIRNAME.
% The file name is guaranteed to be unique.
% FILEID is the handle of the newly opened file.
% FILENAME is the name of the newly opened file. FILENAME
% will consist of 'temp_' followed by a string of numbers repre-
% senting the current date and time to the nearest microsecond.
% If the file cannot be opened (e.g. because your specified
% directory does not exist or you don't have write access
% to it) then an error will occur.
% If only one return argument is used, then the file is created,
% but not opened, and only the filename is returned.
% NOTE: This function differs from Matlab's TEMPNAME function in
% that it is multi-thread safe (i.e. no race-condition), and you
% can specify your own directory instead of /tmp.

% Copyright (C) 2006  Murphy Lab
% Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.
%
% For additional information visit http://murphylab.web.cmu.edu or
% send email to murphy@cmu.edu

% If no input argument, assume /tmp

% Written by Meel Velliste

if( nargin < 1)
    dirname = '/tmp/';
end

% Append / to dirname if it does not already have it
len = length(dirname);
if( dirname(len) ~= '/') dirname(len+1) = '/'; end;

% Try to create a file until successful
while( 1)
    % Generate filename from current time to the nearest microsecond
    c = clock;
    timestr = sprintf('%.4d%.2d%.2d%.2d%.2d%.8d',c(1),c(2),c(3),c(4),c(5),floor(1e6*(c(6))));
    FileName = ['temp_' timestr];
    FullName = [dirname FileName];
    % Create the file if possible
    status = ml_createfile_excl( FullName);
    switch( status)
    case 0,
        if( nargout < 2)
            arg1 = FileName;
        else
            arg1 = fopen( FullName, 'w');
            arg2 = FileName;
        end
        break;
    case 1,
        %Try again
    case 2,
        error(['Could not create temporar file in: ' dirname]);
    otherwise
        error('Unexpected return value from ml_createfile_excl()');
    end;
end;

