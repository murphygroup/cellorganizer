function message = create_comment_header( str )
% CREATE_COMMENT_HEADER Helper function that formats a string to use as a
% header for a comment section

% Ivan E. Cao-Berg
%
% Copyright (C) 2018 Murphy Lab
% Lane Center for Computational Biology
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
% For additional information visit http://murphylab.web.cmu.edu or
% send email to murphy@cmu.edu

if nargin ~=1 | isempty( str ) | ~isa( str, 'char' )
	return
end

str = [ ' ' str ' ' ];
stars(1:ceil((75-length(str))/2))='%';
str = [ stars str stars ];
str = str(1:74);
message = str;
