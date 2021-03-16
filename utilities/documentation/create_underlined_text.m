function message = create_underlined_text( message )
% Ivan E. Cao-Berg
%
% Copyright (C) 2020 Murphy Lab
% Computational Biology Department
% School of Computer Science
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
% For additional information visit http://murphylab.web.cmu.edu/ or
% send email to murphy@cmu.edu

if nargin ~=1 || isempty( message ) || ~isa( message, 'char' )
	return
end

title_length = 75;
if ~isempty( message )
    message = ['  ' message '  '];
    message = [ message; ... 
        repmat( '%', 1, length(message)) ];
end
end%create_underlined_text
