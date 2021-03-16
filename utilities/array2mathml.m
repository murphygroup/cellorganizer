function array2mathml( array, file, tab )
%ARRAY2MATHML Parses a matrix or vector to MathML presentation format.
%
%SYNTAX
%array2mathml( array, file, tab )
%
%COMMENT: This version does not support non-numeric matrixes or vectors

% Author: Ivan E. Cao-Berg (icaoberg@cmu.edu)
% Created: May 8, 2007
% Last Update: March 11, 2008
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

if( nargin == 1 )
    [fileID tabs] = checkInputArguments( array );
    print2screen( array );
elseif( nargin == 3 )
    [fileID tabs] = checkInputArguments( array, tab, file );
    print2file( array, tabs, fileID );
else
    error( 'CellOrganizer: Wrong number of input arguments.' );
end

end%array2mathml

%--------------------------------------------------------------------------
function [fileID tabs] = checkInputArguments( array, tab, file )
% CHECKARGUMENTSIN Checks the input arguments

if( nargin == 1 )
    %check array
    if( ~isnumeric( array ) )
        error( ['CellOrganizer: This method does not support' ...
            'non-numeric arrays'] );
    end
    
    fileID = '';
    tabs = '';
else
    %check array
    if( ~isnumeric( array ) )
        error( ['CellOrganizer: This method does not support' ...
            'non-numeric arrays'] );
    end

    %check tab
    tabs = '';
    if( isempty(tab) )
        return;
    elseif( ischar(tab) )
        tabs = tab;
    elseif( ceil(tab) < 0 )
        error('CellOrganizer: Input argument tab must be nonnegative');
    elseif( tab == 0  )
        tabs = '';
    else
        for i=1:1:tab
            tabs = [ tabs '\t' ];
        end
    end

    %check file
    fileID = file2fileID( file );
end
end%checkInputArguments