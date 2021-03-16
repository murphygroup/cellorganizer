function write( file, name, value, tab, varargin )
%WRITE Simplifies the writing of XML snippets to a file.
%
% Arguments		Description
% ---------     -----------
% file          Filename or file ID
% name          Name of the XML tag
% tab           Number of tabs when writing to file
% value         Value of the XML tag, can be numeric
%               string and/or empty
% varargin      Arguments of the XML tag, these are optional
%	            and strings
%
% Example 1
% ---------
% filename = 'example.xml';
% name = 'mu';
% value = rand(2);
% tab = 0;

% Author: Ivan E. Cao-Berg (icaoberg@cmu.edu)
% Created: May 8, 2007
% Last Update: February 18, 2008
%
% Copyright (C) 2007  Murphy Lab
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
% For additional information visit http://murphylab.web.cmu.edu/software or
% send email to murphy@cmu.edu

if( nargin == 4 )
    [fileID tabs] = checkInputArguments( file, name, value, tab );

    if( isanarray(value) | isaVector(value) )
        fprintf( fileID, [tabs '%s\n'], ['<' name '>'] );
        array2mathml( value, fileID, [tabs '\t'] );
        fprintf( fileID, [tabs '%s\n'], ['</' name '>'] );
    elseif( isnumeric(value) & ( isanarray(value) | isaVector(value) ) )
        fprintf( fileID, [tabs '%s\n'], ['<' name '>' ...
            array2mathml( value, fileID, [tabs '%s\n'] ) '</' name '>' ] );
    elseif( isnumeric(value) )
        fprintf( fileID, [tabs '%s\n'], ['<' name '>' num2str(value) '</' name '>'] );
    elseif( ischar(value) )
        fprintf( fileID, [tabs '%s'], ['<' name '>' ] );
        fprintf( fileID, '%s\n', [ value '</' name '>' ] );
    end
else
    [fileID tabs]= checkInputArguments( file, name, value, tab );

    fprintf( fileID, [tabs '%s'], ['<' name ' '] );
    for( i=1:2:length(varargin) )
        if( i == length(varargin)-1 )
            fprintf( fileID, '%s', [ varargin{i} '="' varargin{i+1} '"'] );
        else
            fprintf( fileID, '%s', [ varargin{i} '="' varargin{i+1} '" '] );
        end
    end

    if( isanArray(value) | isaVector(value) )
        fprintf( fileID, '%s\n', '>' );
        array2mathml( value, fileID, [tabs '\t'] );
        fprintf( fileID, [tabs '%s\n'], ['</' name '>'] );
    elseif( isnumeric(value) & isempty(value) )
        fprintf( fileID, '%s\n', ['>[]</' name '>'] );
    elseif( isnumeric(value) )
        fprintf( fileID, '%s\n', ['>' num2str(value) '</' name '>'] );
    elseif( ischar(value) & isempty(value) )
        fprintf( fileID, '%s\n', '/>' );
    else
        fprintf( fileID, '%s\n', ['>' value '</' name '>'] );
    end
end
end%write

%---------------------------------------------------------------------------
%HELPER FUNCTIONS
function [fileID tabs] = checkInputArguments( file, name, value, tab )
%CHECKARGUMENTSIN

%check file
if( strcmpi( file, '' ) )
    error('CellOrganizer: Input argument file cannot be empty');
elseif( isnumeric( file ) )
    %assume file is a number, thus it may be a Matlab file ID
    if( isempty( fopen( file ) ) )
        error( ['CellOrganizer: Input argument file refers to a file ' ...
            ' ID which points to a nonexisting file'] );
    else
        fileID = file;
    end
elseif( ischar( file ) )
    %assume file is a string, thus it may be a filename
    fileID = fopen( file, 'w' );
else
  error( ['CellOrganizer: Input argument file is either not a valid ' ...
        'file ID or filename'] );
end

%check name
if( ~ischar(name) )
    error('CellOrganizer: Input argument name must be of type string');
end

%check value
if( ~isnumeric(value) && ~ischar(value) && ~isempty(value) )
    error('CellOrganizer: Input argument value must be of type numeric, string or empty');
end

%check tab
if( ischar(tab) )
    tabs = tab;
elseif( ceil(tab) < 0 )
    error('CellOrganizer: Input argument tab must be a nonnegative number');
elseif( isempty(tab) | ceil(tab) == 0 | tab == 0 )
    tabs = '';
else
    tabs = '';
    for( i=1:1:ceil(tab) )
        tabs = [tabs  '\t'];
    end
end
end%checkInputArguments
