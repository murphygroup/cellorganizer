function result = strbegins(str, prefix)
%STRBEGINS Determine if a string begins with a prefix

% Tests
% -----
% assert(strbegins('a', 'a'))
% assert(~strbegins('a', 'b'))
% assert(strbegins('ab', 'a'))
% assert(~strbegins('ab', 'b'))


% Author: Taraz Buck
%
% Copyright (C) 2019-2020 Murphy Lab
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


if ischar(str) && ischar(prefix)
    result = strcmp(str(1:min(length(prefix), end)), prefix);
elseif ischar(str) && iscell(prefix)
    result = cellfun(@(x)strbegins(str, x), prefix);
elseif iscell(str) && ischar(prefix)
    result = cellfun(@(x)strbegins(x, prefix), str);
else
    error('One or both arguments must be string and the other may be cell array containing strings');
end

end

