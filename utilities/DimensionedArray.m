classdef DimensionedArray
    % Array of dimensioned numbers. Also known as a denominate or dimensionful number.
    
    % Author: Taraz Buck
    
    
    properties
        value
        units_manager
    end
    
    
    methods
        
        function obj = DimensionedArray(initial_value, initial_units)
            if nargin < 1
                initial_value = nan;
            elseif nargin < 2
                initial_units = '';
            end
            
            if isa(initial_value, 'DimensionedArray')
                obj.value = initial_value.value;
                obj.units_manager = UnitsManager(initial_value.units_manager);
                return;
            end
            
            obj.value = initial_value;
            obj.units_manager = UnitsManager(initial_units);
            
            obj = obj.simplify();
        end
        
        function disp(obj)
            disp(obj.value);
            disp(obj.units_manager);
        end
        
        function obj = simplify(obj)
            obj = DimensionedArray(obj);
            obj.value = obj.value * obj.units_manager.constant_scale;
            obj.units_manager.constant_scale = 1;
        end
        
        function r = plus(obj1, obj2)
            r = DimensionedArray(obj1);
            if strcmp(class(obj2), 'char')
                obj2 = DimensionedArray(obj2);
            end
            obj2 = obj2.convertTo(r.units_manager);
            r.value = r.value + obj2.value;
            r = r.simplify();
        end
        
        function r = minus(obj1, obj2)
            r = DimensionedArray(obj1);
            if strcmp(class(obj2), 'char')
                obj2 = DimensionedArray(obj2);
            end
            obj2 = obj2.convertTo(r.units_manager);
            r.value = r.value - obj2.value;
            r = r.simplify();
        end
        
        function r = mtimes(obj1, obj2)
            r = DimensionedArray(obj1);
            if strcmp(class(obj2), 'char')
                obj2 = DimensionedArray(obj2);
            end
            r.value = r.value * obj2.value;
            r.units_manager = r.units_manager * obj2.units_manager;
            r = r.simplify();
        end
        
        function r = times(obj1, obj2)
            r = DimensionedArray(obj1);
            if strcmp(class(obj2), 'char')
                obj2 = DimensionedArray(obj2);
            end
            r.value = r.value .* obj2.value;
            r.units_manager = r.units_manager * obj2.units_manager;
            r = r.simplify();
        end
        
        function r = mrdivide(obj1, obj2)
            r = DimensionedArray(obj1);
            if strcmp(class(obj2), 'char')
                obj2 = DimensionedArray(obj2);
            end
            r.value = r.value / obj2.value;
            r.units_manager = r.units_manager / obj2.units_manager;
            r = r.simplify();
        end
        
        function r = rdivide(obj1, obj2)
            r = DimensionedArray(obj1);
            if strcmp(class(obj2), 'char')
                obj2 = DimensionedArray(obj2);
            end
            r.value = r.value ./ obj2.value;
            r.units_manager = r.units_manager / obj2.units_manager;
            r = r.simplify();
        end
        
    end
end