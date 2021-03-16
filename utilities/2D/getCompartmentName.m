function compartmentName = getCompartmentName( file )
%GETCOMPARTMENTNAME Returns the compartment name from a given filename
%assuming the filename may contain the name of the compartment, empty
%otherwise.

% Author: Ivan E. Cao-Berg (icaoberg@cmu.edu)
% Created: May 8, 2007
% Last Updated: March 2, 2008
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
    if( ~exist( file, 'file' ) )
        error( ['CellOrganizer: The file ' file ' does not exist.'] );
    else
        % Get filename from file
        filename = getFileName( file );
        
        if( findstr( filename, 'Lys' ) || findstr( filename, 'LYS' ) ...
                || findstr( file, 'lys' )  )
            compartmentName = 'lysosome';
        elseif( findstr( filename, 'Endo' ) || findstr( filename, 'ENDO' ) ...
                || findstr( filename, 'endo' )  )
            compartmentName = 'endosome';
        elseif( findstr( filename, 'Golgi' ) || findstr( filename, 'GOLGI' ) ...
                || findstr( filename, 'golgi' )  )
            compartmentName = 'golgi';
        else
            compartmentName = '';
        end
    end
else
    error( 'CellOrganizer: Wrong number of input arguments' );
end
end%getCompartmentName