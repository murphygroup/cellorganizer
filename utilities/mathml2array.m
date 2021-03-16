function array = mathml2array( mathml )
%MATHML2ARRAY Parses a MathML vector or array into a Matlab array.

% Author: Ivan E. Cao-Berg (icaoberg@cmu.edu)
% Created: May 8, 2007
% Last Update: March 4, 2008
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
    mathml = checkInputArguments( mathml );

    if( strfind( mathml{1}, 'vector' ) )
        array = vector( mathml );
    elseif( strfind( mathml{1}, 'array' ))
        array = matrix( mathml );
    else
        error( ['CellOrganizer: The information in the cell array is ' ...
            'not a supported MathML code.'] );
    end
end
end%mathml2array

%--------------------------------------------------------------------------
% HELPER METHODS 
function mathml = checkInputArguments( mathml )
% CHECKINPUTARGUMENTS Checks the validity of the arguments of the main function
if( ~isa( mathml, 'cell' ) )
    if( exist( mathml, 'file' ))
        mathml = textread( mathml, '%s', 'delimiter', '\n' );
    else
        error( 'CellOrganizer: Input must me a cell array or an XML file.' );
    end
end
end%checkInputArguments