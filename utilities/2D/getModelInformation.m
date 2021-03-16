function information = getModelInformation( data )
%GETDOCUMENTATION Returns the documentation of valid SLML model

% Author: Ivan E. Cao-Berg (icaoberg@cmu.edu)
% Created: May 29, 2007
%
% Copyright (C) 2008  Murphy Lab
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

% Analyze documentation
delimiters = [];
for i=1:1:length( data )
    if( strfind( data{i}, 'slml' ) )
        delimiters = [ delimiters, i ];
    end
end

if( isempty( delimiters ) )
    % No documentation
    documentation = struct([]);
else
    %Preallocate in memory
    documentation.temp = 'temp';

    % In a valid SLML model you only expect a single documentation tag thus
    % expecting to get only two delimiters is valid
    for i=delimiters(1)+1:1:delimiters(2)-1
        quotations = findstr( '"', data{i} );
        documentation = setfield( documentation, ...
            data{i}(quotations(1)+1:quotations(2)-1), ...
            data{i}(quotations(3)+1:quotations(4)-1));
    end
end

% Remove temporary field
documentation = rmfield( documentation, 'temp' );
end%getModelInformation