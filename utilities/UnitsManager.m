classdef UnitsManager
% Manage units or dimensions through arithmetic operations.
% 
% %Example 1: Create a UnitsManager from a string
% UnitsManager('10 uM')
% 
% %Example 2: Multiplying two UnitsManagers
% UnitsManager('10') * UnitsManager('uM')

% Tests
% -----
% assert(UnitsManager('2 * km.h-1') == UnitsManager('5.5555555555555558e-1 * m.s-1'))
% assert(UnitsManager('2 * km.h-1') * UnitsManager('9 * s') == UnitsManager('5 * m'))
% assert(UnitsManager('2 m')^2 == UnitsManager('4 * m2'))
% assert(UnitsManager('1 m') * UnitsManager('1 mN') == UnitsManager('1e-6 * kN.m'))

% Author: Taraz Buck
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
    
    
    properties(Constant, Hidden)
        prefixes_full_array = {'yotta', 'zetta', 'exa', 'peta', 'tera', 'giga', 'mega', 'kilo', 'hecto', 'deca', 'deci', 'centi', 'milli', 'micro', 'nano', 'pico', 'femto', 'atto', 'zepto', 'yocto'};
        prefixes_abbreviated_array = {'Y', 'Z', 'E', 'P', 'T', 'G', 'M', 'k', 'h', 'da', 'd', 'c', 'm', 'u', 'n', 'p', 'f', 'a', 'z', 'y'};
        prefixes_full_scales_array = [1e24, 1e21, 1e18, 1e15, 1e12, 1e9, 1e6, 1e3, 1e2, 1e1, 1e-1, 1e-2, 1e-3, 1e-6, 1e-9, 1e-12, 1e-15, 1e-18, 1e-21, 1e-24];
        prefixes_array = [UnitsManager.prefixes_full_array, UnitsManager.prefixes_abbreviated_array];
        prefixes_abbreviated_to_full_map = containers.Map(UnitsManager.prefixes_abbreviated_array, UnitsManager.prefixes_full_array);
        prefixes_scales = containers.Map(UnitsManager.prefixes_array, [UnitsManager.prefixes_full_scales_array, UnitsManager.prefixes_full_scales_array]);
        
        % Use grams instead of kilograms to simplify implementation. Probably does not make a difference in practical use.
        base_units_full_array = {'metre', 'gram', 'second', 'ampere', 'kelvin', 'mole', 'candela', 'molecule'};
        base_units_abbreviated_to_full_map = containers.Map({'m', 'meter', 'g', 'sec', 's', 'amp', 'A', 'K', 'mol', 'cd'}, {'metre', 'metre', 'gram', 'second', 'second', 'ampere', 'ampere', 'kelvin', 'mole', 'candela'});
        base_units_array = [UnitsManager.base_units_full_array, sort(UnitsManager.base_units_abbreviated_to_full_map.keys())];
        
        derived_units_to_base_units_map = containers.Map({'liter', 'L', 'M', 'N', 'V', 'h'}, {'dm3', 'dm3', 'mol.dm-3', 'kg.m.s-2', 'kg.m2.s-3.A-1', '3600 * s'});
        
        all_units_array = [UnitsManager.base_units_array, sort(UnitsManager.derived_units_to_base_units_map.keys())];
        
        prefix_pattern = ['(', strjoin(UnitsManager.prefixes_array, '|'), ')'];
        unit_pattern = ['(', strjoin(UnitsManager.all_units_array, '|'), ')'];
        exponent_pattern = ['(-?[0-9]+)'];
        complete_unit_pattern = ['(?<prefix>(', UnitsManager.prefix_pattern, ')?)', '(?<base_unit>(', UnitsManager.unit_pattern, '))', '(?<exponent>(', UnitsManager.exponent_pattern, ')?)'];
        unit_delimiter = '.';
        
        simple_float_pattern = ['[-+]?([0-9]*\.?[0-9]+|[0-9]+\.?[0-9]*)([Ee][-+]?[0-9]+)?'];
        constant_delimiters = {' ', '*'};
    end
    
    
    properties
        unit_exponents
        constant_scale
    end
    
    
    methods(Static)
        
        function [prefix, base_unit, exponent] = splitUnit(unit_str)
            named_tokens = regexp(unit_str, ['^', UnitsManager.complete_unit_pattern, '$'], 'names');
            prefix = named_tokens.prefix;
            base_unit = named_tokens.base_unit;
            exponent = named_tokens.exponent;
        end
        
        function [base_unit, exponent, scale] = parseUnit(unit_str)
            [prefix, base_unit, exponent] = UnitsManager.splitUnit(unit_str);
            if isempty(exponent)
                exponent = 1;
            else
                exponent = str2double(exponent);
            end
            if isempty(prefix)
                scale = 1;
            else
                scale = UnitsManager.prefixes_scales(prefix)^exponent;
            end
        end
        
        function [unit_strs] = split(units_str)
            unit_strs = strsplit(units_str, UnitsManager.unit_delimiter);
        end
        
        function r = baseUnit(unit)
            obj.units_info = containers.Map('KeyType', 'char', 'ValueType', 'any');
            obj.constant_scale = 1;
        end
        
        function [obj1, obj2] = processArithmeticArguments(obj1, obj2)
            obj1 = UnitsManager(obj1);
            obj2 = UnitsManager(obj2);
        end
        
    end
    
    
    methods
        
        function obj = UnitsManager(units_str)
            if nargin < 1
                units_str = 1;
            end
            
            if isa(units_str, 'UnitsManager')
                if length(units_str.unit_exponents) > 0
                    obj.unit_exponents = containers.Map(units_str.unit_exponents.keys(), units_str.unit_exponents.values());
                else
                    obj.unit_exponents = containers.Map('KeyType', 'char', 'ValueType', 'double');
                end
                obj.constant_scale = units_str.constant_scale;
                return;
            end
            
            if isnumeric(units_str)
                obj.unit_exponents = containers.Map('KeyType', 'char', 'ValueType', 'double');
                obj.constant_scale = double(units_str);
                return;
            end
            
            if ~ischar(units_str)
                error('initial_value must be numeric, string, or UnitsManager');
            end
            
            obj.unit_exponents = containers.Map('KeyType', 'char', 'ValueType', 'double');
            obj.constant_scale = 1;
            
            units_str = strtrim(units_str);
            
            % Process optional constant
            units_str_parts = strsplit(units_str, UnitsManager.constant_delimiters);
            units_str_parts = strtrim(units_str_parts);
            %{
            units_str_float_part = regexp(units_str_parts{1}, ['^', UnitsManager.simple_float_pattern, '$'], 'match');
            if length(units_str_parts) > 1
                constant_str = units_str_parts{1};
                units_str = units_str_parts{2};
                constant_str = strtrim(constant_str);
                units_str = strtrim(units_str);
                obj.constant_scale = obj.constant_scale * str2double(constant_str);
            elseif ~isempty(units_str_float_part)
                obj.constant_scale = obj.constant_scale * str2double(units_str_float_part);
                units_str = '';
            end
            %}
            if length(units_str_parts) == 0 || length(units_str_parts) > 2
                error('initial_value string must contain only a number, units, or both');
            end
            units_str_constant_part = '';
            if strcmpi(units_str_parts{1}, 'nan') || ~isnan(str2double(units_str_parts{1}))
                units_str_constant_part = units_str_parts{1};
                units_str_parts = units_str_parts(2);
                obj.constant_scale = obj.constant_scale * str2double(units_str_constant_part);
            end
            units_str_units_part = '';
            if length(units_str_parts) == 1
                units_str_units_part = units_str_parts{1};
                units_str = units_str_units_part;
            end
            
            % Process units
            unit_strs = UnitsManager.split(units_str);
            if isempty(units_str)
                unit_strs = {};
            end
            for i = 1:length(unit_strs)
                unit_str = unit_strs{i};
                [base_unit, exponent, scale] = UnitsManager.parseUnit(unit_str);
                obj.constant_scale = obj.constant_scale * scale;
                if UnitsManager.derived_units_to_base_units_map.isKey(base_unit)
                    temp = UnitsManager(UnitsManager.derived_units_to_base_units_map(base_unit))^exponent;
                    obj.constant_scale = obj.constant_scale * temp.constant_scale;
                    temp_unit_exponents = temp.unit_exponents;
                    temp_base_units = temp.unit_exponents.keys();
                    for j = 1:length(temp_unit_exponents)
                        temp_base_unit = temp_base_units{j};
                        % temp_exponent = temp_unit_exponents(temp_base_unit) * exponent;
                        temp_exponent = temp_unit_exponents(temp_base_unit);
                        if ~obj.unit_exponents.isKey(temp_base_unit)
                            obj.unit_exponents(temp_base_unit) = 0;
                        end
                        obj.unit_exponents(temp_base_unit) = obj.unit_exponents(temp_base_unit) + temp_exponent;
                    end
                else
                    if ~obj.unit_exponents.isKey(base_unit)
                        obj.unit_exponents(base_unit) = 0;
                    end
                    obj.unit_exponents(base_unit) = obj.unit_exponents(base_unit) + exponent;
                end
            end
            
            obj = obj.simplify();
        end
        
        function r = isDimensionless(obj)
            r = length(obj.unit_exponents) == 0;
        end
        
        function r = toString(obj, with_scale, format_variant)
            if nargin < 2
                with_scale = true;
            end
            if nargin < 3
                format_variant = '';
            end
            obj = obj.simplify();
            
            r = '';
            units_keys = sort(obj.unit_exponents.keys());
            
            if with_scale
                r = [r, double2str(obj.constant_scale)];
                if length(units_keys) > 0
                    if UnitsManager.constant_delimiters{1} == ' '
                        r = [r, ' '];
                    else
                        r = [r, ' ', UnitsManager.constant_delimiters{1}, ' '];
                    end
                end
            end
            
            for i = 1:length(units_keys)
                if i > 1
                    r = [r, '.'];
                end
                units_key = units_keys{i};
                if strcmp(format_variant, 'VCML')
                    if strcmp(units_key, 'molecule')
                        units_key = 'molecules';
                    end
                end
                r = [r, units_key];
                exponent = obj.unit_exponents(units_keys{i});
                if exponent ~= 1
                    r = [r, num2str(exponent)];
                end
            end
            if strcmp(format_variant, 'VCML')
                if length(units_keys) == 0
                    r = '1';
                end
            end
        end
        
        function r = char(obj)
            r = obj.toString();
        end
        
        function disp(obj)
            disp(['    ', obj.toString()]);
            fprintf('\n');
        end
        
        function r = double(obj)
            obj = obj.simplify();
            if length(obj.unit_exponents) > 0
                error('UnitsManager has nonzero unit exponents');
            end
            r = obj.constant_scale;
        end
        
        function obj = simplify(obj)
            obj = UnitsManager(obj);
            
            obj_unit_exponents_keys = obj.unit_exponents.keys();
            for i = 1:length(obj_unit_exponents_keys)
                obj_unit = obj_unit_exponents_keys{i};
                exponent = obj.unit_exponents(obj_unit);
                if exponent == 0
                    obj.unit_exponents.remove(obj_unit);
                end
            end
        end
        
        function r = eq(obj1, obj2, ignore_scale)
            % Check if units match (scales can be different if ignore_scale is true)
            if nargin < 3
                ignore_scale = false;
            end
            [obj1, obj2] = UnitsManager.processArithmeticArguments(obj1, obj2);
            
            obj2_unit_exponents_keys = obj2.unit_exponents.keys();
            if length(obj1.unit_exponents.keys()) ~= length(obj2_unit_exponents_keys)
                r = false;
                return;
            end
            for i = 1:length(obj2_unit_exponents_keys)
                key = obj2_unit_exponents_keys{i};
                if ~obj1.unit_exponents.isKey(key)
                    r = false;
                    return;
                end
                if obj1.unit_exponents(key) ~= obj2.unit_exponents(key)
                    r = false;
                    return;
                end
            end
            if ~ignore_scale && obj1.constant_scale ~= obj2.constant_scale
                r = false;
                return;
            end
            r = true;
        end
        
        function r = ne(obj1, obj2)
            r = ~eq(obj1, obj2);
        end
        
        function r = mtimes(obj1, obj2)
            [r, obj2] = UnitsManager.processArithmeticArguments(obj1, obj2);
            
            obj2_unit_exponents_keys = obj2.unit_exponents.keys();
            for i = 1:length(obj2_unit_exponents_keys)
                key = obj2_unit_exponents_keys{i};
                if ~r.unit_exponents.isKey(key)
                    r.unit_exponents(key) = 0;
                end
                r.unit_exponents(key) = r.unit_exponents(key) + obj2.unit_exponents(key);
            end
            r.constant_scale = r.constant_scale * obj2.constant_scale;
            
            r = r.simplify();
        end
        
        function r = times(obj1, obj2)
            r = mtimes(obj1, obj2);
        end
        
        function r = mrdivide(obj1, obj2)
            [r, obj2] = UnitsManager.processArithmeticArguments(obj1, obj2);
            
            obj2_unit_exponents_keys = obj2.unit_exponents.keys();
            for i = 1:length(obj2_unit_exponents_keys)
                key = obj2_unit_exponents_keys{i};
                if ~r.unit_exponents.isKey(key)
                    r.unit_exponents(key) = 0;
                end
                r.unit_exponents(key) = r.unit_exponents(key) - obj2.unit_exponents(key);
            end
            r.constant_scale = r.constant_scale / obj2.constant_scale;
            
            r = r.simplify();
        end
        
        function r = rdivide(obj1, obj2)
            r = mrdivide(obj1, obj2);
        end
        
        function r = power(obj1, val2)
            if isinteger(val2)
                val2 = double(val2);
            end
            if ~isfloat(val2)
                error('Second argument must be integer or float representing compatible fraction');
            end
            [val2_num, val2_denom] = rat(val2);
            
            [r, ~] = UnitsManager.processArithmeticArguments(obj1, val2);
            
            r_unit_exponents_keys = r.unit_exponents.keys();
            for i = 1:length(r_unit_exponents_keys)
                key = r_unit_exponents_keys{i};
                r.unit_exponents(key) = r.unit_exponents(key) * val2_num / val2_denom;
                if abs(r.unit_exponents(key) - round(r.unit_exponents(key))) > 1e-5
                    error('Second argument must be integer or float representing compatible fraction');
                end
                r.unit_exponents(key) = round(r.unit_exponents(key));
            end
            r.constant_scale = r.constant_scale^val2;
            
            r = r.simplify();
        end
        
        function r = mpower(obj1, val2)
            r = power(obj1, val2);
        end
        
    end
    
end
