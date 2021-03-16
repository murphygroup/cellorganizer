function print2file( array, tabs, fileID )
%PRINT2FILE Helper function that prints the contents of an array on a
%file using MathML presentation format.

% Author: Ivan E. Cao-Berg (icaoberg@cmu.edu)
% Created: May 8, 2007
% Last Update: March 19, 2008
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

%to avoid truncation
format long

% print an array
fprintf( fileID, [tabs '%s\n'], ['<array dimensions="(' num2str(size(array,1)) ',' num2str(size(array,2)) ')">'] );
for i=1:1:size(array,1)
    for j=1:1:size(array,2)
        if array(i,j) ~= 0
            fprintf( fileID, [tabs '\t%s\n'], ['<number index="(' num2str(i) ','  num2str(j) ')">' num2str(array(i,j)) '</number>'] );
        end
    end
end
fprintf( fileID, [tabs '%s\n'], '</array>' );
end%print2file
