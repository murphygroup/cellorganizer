function files = getFiles( directory, option )
%GETFILES Returns a cell array containing all the filenames and/or paths in
%the specified directory.
%
%Options
%filename            Returns a cell array containing only the filenames
%fullpath            Returns a cell array containing the filenamew and paths
%
%Example 1
%> getFiles; %Returns all the files in the current folder
%
%Example 2
%> getFiles( pwd ); %The same as Example 1
%
%Example 3
%> directory = '~/data';
%> getFiles( directory, 'fullpath' ); %Returns the filenames and path

% Author: Ivan E. Cao-Berg (icaoberg@cmu.edu)
% Created: September 17, 2007
% Last Update: March 2, 2008
%
% Copyright (C) 2008 Center for Bioimage Informatics/Murphy Lab
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

files = [];
%getFiles( pwd )
if( nargin == 0 )
    directory = pwd;
    directoryInfo = dir( directory );
    for i=3:1:length(directoryInfo)
        if( isaFile( [directory filesep directoryInfo(i).name] ) )
            files{length(files)+1} = directoryInfo(i).name; %#ok<AGROW>
        end
    end
% getFiles( directory )
elseif( nargin == 1 )
    if( ~exist( directory, 'dir' ) )
        warning( 'CellOrganizer: Directory nonexistent or invalid' ); %#ok<WNTAG>
        return
    else
        directoryInfo = dir( directory );
        files = [];
        for i=3:1:length(directoryInfo)
            if( isaFile( [directory filesep directoryInfo(i).name] ))
                files{length(files)+1} = [ directory filesep directoryInfo(i).name ]; %#ok<AGROW>
            end
        end
    end
    % getFiles( directory, option )
elseif( nargin == 2 )
    if( ~isdir( directory ) )
        error( 'CellOrganizer: Directory nonexistent or invalid' ); %#ok<WNTAG>
        % getFiles( directory, 'filename' )
    elseif( strcmpi( option, 'filename' ) )
        directoryInfo = dir( directory );
        for i=3:1:length(directoryInfo)
            if( isaFile( [directory filesep directoryInfo(i).name] ) )
                files{length(files)+1} = directoryInfo(i).name; %#ok<AGROW>
            end
        end
        % getFiles( directory, 'fullpath' )
    elseif( strcmpi( option, 'fullpath' ) )
        directoryInfo = dir( directory );
        for i=3:1:length(directoryInfo)
            if( isaFile( [directory filesep directoryInfo(i).name] ) )
                files{length(files)+1} = [ directory filesep directoryInfo(i).name ]; %#ok<AGROW>
            end
        end
        % getFiles( directory, 'anything_else' )
    else
        error('CellOrganizer: Unsupported option for this function' );
    end
else
    error( 'CellOrganizer: Wrong number of input arguments.' );
end
 end%getFiles