function result = strends(str, suffix)
%STRENDS Determine if a string ends with a suffix

% Tests
% -----
% assert(strends('a', 'a'))
% assert(~strends('a', 'b'))
% assert(strends('ba', 'a'))
% assert(~strends('ba', 'b'))


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


if ischar(str) && ischar(suffix)
    result = strcmp(str(max(end-length(suffix)+1, 1):end), suffix);
elseif ischar(str) && iscell(suffix)
    result = cellfun(@(x)strends(str, x), suffix);
elseif iscell(str) && ischar(suffix)
    result = cellfun(@(x)strends(x, suffix), str);
else
    error('One or both arguments must be string and the other may be cell array containing strings');
end

end

