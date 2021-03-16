function documentation = getDocumentation( filename )
%GETDOCUMENTATION Returns the documentation of an SLML instance as a Matlab
%structure.

% Author: Ivan E. Cao-Berg (icaoberg@cmu.edu)
% Created: May 29, 2007
% Last Update: March 8, 2008
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

if( nargin ~= 1 )
    error( 'CellOrganizer: Wrong number of input arguments' );
else
    if( ~isaFile( filename ) )
        error( 'CellOrganizer: Invalid input argument' );
    else
        %get data and save as a cell array
        data = file2cell( filename );
        documentation = struct([]);

        %get the documentation section delimiters
        delimiters = getDocumentationDelimiters( data );

        if ~isempty(delimiters)
            %simplify data
            data = data(delimiters(1):delimiters(2));

            %parse data as a Matlab structure
            for i=2:1:length(data)-1
                documentation(i-1).name = getName( data{i} );
                documentation(i-1).value = getValue( data{i} );
            end
        else
            %no documentation found in the SLML instance
            return
        end
    end
end
end%getDocumentation

%--------------------------------------------------------------------------
function name = getName( snippet )
%GETNAME Helper function that searches and returns the name of the current
%XML tag
delimiters = snippet( string, '"' );
name = snippet(delimiters(1)+1:delimiters(2)-1);
end%getName

function value = getValue( snippet )
%GETVALUE Helper function that searches and returns the value of the current
%XML tag
delimiters = findstr( snippet, '"' );
value = snippet(delimiters(3)+1:delimiters(4)-1);
end%getValue