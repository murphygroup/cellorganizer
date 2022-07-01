function result = isfieldr(param, varargin)
%ISFIELDR Determine if subfield of struct exists

% Tests
% -----
% a = struct(); a.b.c.d = 1; assert(isfieldr(a, 'b.c.d')); assert(isfieldr(a, 'b', 'c.d')); assert(~isfieldr(a, 'b.c.d.e'))


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


if nargin < 2
    error('At least 2 arguments are required')
end

debug_mode = false;
% debug_mode = true;
if debug_mode; fprintf('debug_mode = ''%d''\n', debug_mode); end

fieldpath_components = cell(0, 1);
for i = 1:length(varargin)
    fieldpath_components2 = split(varargin{i}, '.');
    fieldpath_components = [fieldpath_components; fieldpath_components2];
end
if debug_mode; fieldpath_components, end
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
result = true;
