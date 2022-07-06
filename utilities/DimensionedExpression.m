classdef DimensionedExpression
% Expression with units. Also known as a denominate or dimensionful number.

% Tests
% -----
% assert(DimensionedExpression(2, 'km.h-1') * DimensionedExpression(5, 'h') == DimensionedExpression(10, 'km'))

% Author: Taraz Buck
    
    
    properties
        value
        units_manager
    end
    
    
    methods(Static)
        
        function [obj1, obj2] = processArithmeticArguments(obj1, obj2)
            obj1 = DimensionedExpression(obj1);
            obj2 = DimensionedExpression(obj2);
        end
        
    end
    
    
    methods
        
        function obj = DimensionedExpression(initial_value, initial_units)
            if nargin < 1
                initial_value = '0';
            end
            if nargin < 2
                initial_units = '';
            end
            
            if isa(initial_value, 'DimensionedExpression')
                temp = initial_value;
                initial_value = temp.value;
                initial_units = UnitsManager(temp.units_manager) * UnitsManager(initial_units);
            elseif isa(initial_value, 'UnitsManager')
                temp = initial_value;
                initial_value = '1';
                initial_units = UnitsManager(temp) * UnitsManager(initial_units);
            elseif isnumeric(initial_value)
                initial_value = double2str(initial_value);
                initial_units = UnitsManager(initial_units);
            end
            
            if ~ischar(initial_value)
                error('initial_value must be numeric, string, DimensionedExpression, or UnitsManager');
            end
            
            obj.value = initial_value;
            obj.units_manager = UnitsManager(initial_units);
            
            % fprintf('*** obj.value = %s\n', obj.value);
            obj = obj.simplify(true);
            % fprintf('*** obj.value after simplify = %s\n', obj.value);
        end
        
        function obj = set.value(obj, value)
            if isnumeric(value)
                obj.value = double2str(value);
            elseif ischar(value)
                obj.value = value;
            else
                error('value must be numeric or string');
            end
        end
        
        function obj = set.units_manager(obj, units_manager)
            if isa(units_manager, 'UnitsManager') || ischar(units_manager)
                units_manager = UnitsManager(units_manager);
                if units_manager.constant_scale ~= 1
                    obj.value = sprintf('(%s) * (%s)', obj.value, double2str(units_manager.constant_scale));
                    units_manager.constant_scale = 1;
                end
                obj.units_manager = units_manager;
            else
                error('value must be UnitsManager or string');
            end
        end
        
        function r = isDimensionless(obj)
            r = obj.units_manager.isDimensionless();
        end
        
        function r = toString(obj)
            if numel(obj) == 1
                obj = obj.simplify();
                r = obj.value;
                if length(obj.units_manager.unit_exponents) > 0 || obj.units_manager.constant_scale ~= 1
                    r = [r, ' * ', obj.units_manager.toString()];
                end
            else
                r = sprintf('%dx%d DimensionedExpression array', size(obj));
            end
        end
        
        function r = char(obj)
            r = obj.toString();
        end
        
        function disp(obj)
            disp(['    ', obj.toString()]);
            fprintf('\n');
        end
        
        function obj = simplify(obj, in_place)
            if nargin < 2
                in_place = false;
            end
            % fprintf('*** in_place = %i\n', in_place);
            if ~in_place
                % fprintf('*** Replacing\n');
                obj = DimensionedExpression(obj);
            end
            if ~obj.hasVariables()
                % fprintf('*** Evaluating\n');
                obj.value = mathEval(obj.value, [], struct('check_variables', false, 'return_type', 'double'));
            end
        end
        
        function r = eq(obj1, obj2)
            % Check if expression and units match
            [obj1, obj2] = DimensionedExpression.processArithmeticArguments(obj1, obj2);
            
            if ~strcmp(obj1.value, obj2.value)
                r = false;
                return;
            end
            if obj1.units_manager ~= obj2.units_manager
                r = false;
                return;
            end
            r = true;
        end
        
        function r = ne(obj1, obj2)
            r = ~eq(obj1, obj2);
        end
        
        function obj = replaceVariable(obj, old_variable, new_expression)
            obj = DimensionedExpression(obj);
            obj.value = regexprep(obj.value, ['\<', old_variable, '\>'], new_expression);
        end
        
        function value_tokens = getVariables(obj)
            % obj = DimensionedExpression(obj);
            value_tokens = mathTokenize(obj.value, [], struct('check_variables', false));
            value_tokens = unique({value_tokens(strcmp({value_tokens.type}, 'variable')).token});
        end
        
        function r = hasVariables(obj)
            % obj = DimensionedExpression(obj);
            r = length(obj.getVariables()) > 0;
        end
        
        function r = plus(obj1, obj2)
            [r, obj2] = DimensionedExpression.processArithmeticArguments(obj1, obj2);
            if str2double(obj2.value) == 0
                r = r;
            elseif str2double(r.value) == 0
                r = obj2;
            else
                if ~eq(obj1.units_manager, obj2.units_manager, true)
                    error('Units do not match');
                end
                r.units_manager.constant_scale = obj2.units_manager.constant_scale;
                if obj1.units_manager.constant_scale == obj2.units_manager.constant_scale
                    r.value = sprintf('(%s) + (%s)', r.value, obj2.value);
                else
                    r.value = sprintf('((%s) * (%f / %f)) + (%s)', r.value, obj1.units_manager.constant_scale, obj2.units_manager.constant_scale, obj2.value);
                end
            end
            r = r.simplify();
        end
        
        function r = minus(obj1, obj2)
            [r, obj2] = DimensionedExpression.processArithmeticArguments(obj1, obj2);
            if str2double(obj2.value) == 0
                r = r;
            elseif str2double(r.value) == 0
                r = -1 * obj2;
            else
                if ~eq(obj1.units_manager, obj2.units_manager, true)
                    error('Units do not match');
                end
                r.units_manager.constant_scale = obj2.units_manager.constant_scale;
                if obj1.units_manager.constant_scale == obj2.units_manager.constant_scale
                    r.value = sprintf('(%s) - (%s)', r.value, obj2.value);
                else
                    r.value = sprintf('((%s) * (%f / %f)) - (%s)', r.value, obj1.units_manager.constant_scale, obj2.units_manager.constant_scale, obj2.value);
                end
            end
            r = r.simplify();
        end
        
        function r = mtimes(obj1, obj2)
            [r, obj2] = DimensionedExpression.processArithmeticArguments(obj1, obj2);
            r.units_manager = r.units_manager * obj2.units_manager;
            if str2double(r.value) == 0 && str2double(obj2.value) == 0
                r.value = '0';
            elseif str2double(r.value) == 1
                r.value = obj2.value;
            elseif str2double(obj2.value) == 1
                r.value = r.value;
            else
                r.value = sprintf('(%s) * (%s)', r.value, obj2.value);
            end
            r = r.simplify();
        end
        
        function r = times(obj1, obj2)
            r = mtimes(obj1, obj2);
        end
        
        function r = mrdivide(obj1, obj2)
            [r, obj2] = DimensionedExpression.processArithmeticArguments(obj1, obj2);
            r.units_manager = r.units_manager / obj2.units_manager;
            if str2double(r.value) == 0 && str2double(obj2.value) == 0
                r.value = 'nan';
            elseif str2double(obj2.value) == 1
                r.value = r.value;
            else
                r.value = sprintf('(%s) / (%s)', r.value, obj2.value);
            end
            r = r.simplify();
        end
        
        function r = rdivide(obj1, obj2)
            r = mrdivide(obj1, obj2);
        end
        
        function r = power(obj1, val2)
            if isinteger(val2)
                val2 = double(val2);
            end
            if isa(val2, 'DimensionedExpression') && length(val2.getVariables()) == 0
                val2 = str2double(val2.value);
            end
            
            [r, ~] = DimensionedExpression.processArithmeticArguments(obj1, val2);
            
            r.units_manager = r.units_manager ^ val2;
            if str2double(r.value) == 0 && val2 == 0
                r.value = 'nan';
            elseif (str2double(r.value) > 0 || str2double(r.value) < 0) && val2 == 0
                r.value = '1';
            elseif str2double(r.value) == 0 && (str2double(r.value) > 0 || str2double(r.value) < 0)
                r.value = '0';
            elseif str2double(r.value) == 1 && (str2double(r.value) > 0 || str2double(r.value) < 0)
                r.value = '1';
            else
                [val2_num, val2_denom] = rat(val2);
                if val2_denom > 1
                    r.value = sprintf('(%s) ^ (%d/%d)', r.value, val2_num, val2_denom);
                else
                    r.value = sprintf('(%s) ^ %d', r.value, val2_num);
                end
            end
            r = r.simplify();
        end
        
        function r = mpower(obj1, val2)
            r = power(obj1, val2);
        end
        
        function r = prod(objs)
            if length(objs) == 0
                r = DimensionedExpression('0');
                return;
            end
            r = DimensionedExpression(objs(1));
            for i = 2:length(objs)
                r = r * objs(i);
            end
        end
        
    end
end
