function result = unitsPower(units_string, exponent)
%UNITSPOWER Raises units string units_string to power exponent
%
% Inputs
% ------
% units_string = string representing units with powers, joined with "."
% exponent     = integer
%
% Outputs
% -------
% result = boolean flag indicating success
%
% Notes
% -----
% * 


% Authors: Taraz Buck
%
% Copyright (C) 2019 Murphy Lab
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

% 2019-07-09 Taraz Buck - Created


units_string_tokens = regexp(units_string, '([A-Za-z]+)(-?[0-9]+)?', 'tokens');
result = '';
for i = 1:length(units_string_tokens)
    unit_string = units_string_tokens{i}{1};
    unit_exponent = str2num(units_string_tokens{i}{2});
    if isempty(unit_exponent)
        unit_exponent = 1;
    end
    unit_exponent = unit_exponent * exponent;
    if unit_exponent == 0
        continue;
    end
    if length(result) > 0
        result = [result, '.'];
    end
    result = [result, unit_string];
    if unit_exponent ~= 1
        result = [result, num2str(unit_exponent)];
    end
end

