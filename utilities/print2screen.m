function print2screen( array )
%PRINT2SCREEN Helper function that prints the contents of an array on the
%screen using MathML presentation format.

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

arraySize = size( array );

% print a row vector
if( arraySize(1) == 1 )
    disp( '<vector>' );
    for j=1:1:arraySize(2)
        disp( [ '<cn>' num2str(array(1,j)) '</cn>' ] );
    end
    disp( '</vector>' );
    % print a column vector
elseif( arraySize(2) == 1 )
    disp( '<vector>' );
    for i=1:1:arraySize(1)
        disp( [ '<cn>' num2str(array(i,1)) '</cn>' ] );
    end
    disp( '</vector>' );
    % print an array
else
    disp( '<array>' );
    for i=1:1:arraySize(1)
        disp( '<arrayrow>' );
        for j=1:1:arraySize(2)
            disp( [ '<cn>' num2str(array(i,j)) '</cn>' ] );
        end
        disp( '</arrayrow>' );
    end
    disp( '</array>' );
end
end%print2screen