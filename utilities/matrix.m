function array = matrix( mathml )
%MATRIX Helper function that handles the reading of MathML matrices in
%presentation format.

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

noRows = 0;
noColumns = 0;

% Calculate the number of rows
for i=2:1:length(mathml)-1
    if( findstr( mathml{i}, 'arrayrow' ) )
        noRows = noRows + 1;
    end
end
noRows = noRows/2;

% Calculate the number of columns
for i=2:1:length(mathml)-1
    if( findstr( mathml{i}, 'cn' ) )
        noColumns = noColumns + 1;
    end
end
noColumns = noColumns / noRows;

% Create an empty matrix of fixed size
array = zeros( noRows, noColumns );

% Fill the empty matrix
temp = [];
for i=2:1:length(mathml)-1
    if( findstr( mathml{i}, 'cn' ));
        delimiter = findstr( mathml{i}, '<' );
        temp = [temp, str2double(mathml{i}(5:delimiter(2)-1))];
    end
end

counter = 1;
for i=1:1:noRows
    for j=1:1:noColumns
        array(i,j)=temp(counter);
        counter = counter + 1;
    end
end
end%matrix