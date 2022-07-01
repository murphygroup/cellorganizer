function result = fieldisequal(param, fieldpath, fieldvalue, casesensitive)
%FIELDISEQUAL Determine if subfield of struct exists and ISEQUAL to value

% Tests
% -----
% a = struct(); a.b.c.d = 1; assert(fieldisequal(a, 'b.c.d', 1)); assert(~fieldisequal(a, 'b.c.d', '1')); assert(~fieldisequal(a, 'b.c.d.e', 1))
% a = struct(); a.b.c.d = 'Asdf'; assert(~fieldisequal(a, 'b.c.d', 'asdf')); assert(fieldisequal(a, 'b.c.d', 'asdf', false))
% a = struct(); a.b.c.d = {'A', 'b'}; assert(~fieldisequal(a, 'b.c.d', 'A')); assert(~fieldisequal(a, 'b.c.d', {'a', 'b'})); assert(fieldisequal(a, 'b.c.d', {'A', 'b'}))


% Author: Taraz Buck
%
% Copyright (C) 2022 Murphy Lab
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


if nargin < 3
    error('Exactly 3 arguments are required')
end

debug_mode = false;
% debug_mode = true;
if debug_mode; fprintf('debug_mode = ''%d''\n', debug_mode); end

if nargin < 4
    casesensitive = true;
end
if debug_mode; fprintf('casesensitive = ''%d''\n', casesensitive); end

fieldpath_components = split(fieldpath, '.');
param2 = param;
for i = 1:length(fieldpath_components)
    fieldpath_component = fieldpath_components{i};
    if debug_mode; fprintf('fieldpath_component = ''%s''\n', fieldpath_component); end
    if ~isfield(param2, fieldpath_component)
        result = false;
        return;
    end
    param2 = param2.(fieldpath_component);
end

function result = isequal_sensitive(a, b)
    if ~casesensitive
        if ischar(a) && ischar(b)
            result = strcmpi(a, b);
            return;
        elseif iscellstr(a) && iscellstr(b) && isequal(size(a), size(b))
            result = all(strcmpi(a, b));
            return;
        end
    end
    result = isequal(a, b);
end

% if debug_mode; fprintf('param2 = ''%s''\n', param2); end
if debug_mode; param2, end
if ~isequal_sensitive(param2, fieldvalue)
    result = false;
    return;
end

result = true;

end
