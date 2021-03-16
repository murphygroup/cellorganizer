function [options] = process_options_structure_fast(default_options, options)
% Return a structure the same as the structure options but with fields added from the structure default_options if they are not in options. Produces an error if any string in the cell array of strings required_options is not a field name in options. The second and third arguments are optional. Applies process_options_structure to fields of default_options that are structures if process_substructures = true.
% 
% xruan 09/18/2015

    if nargin == 1
        options = default_options; 
        return
    end

    options_1 = default_options;
    options_names = fieldnames(options);

    for index = 1:length(options_names)
        current_option = options_names{index};
        options_1.(current_option) = options.(current_option);
    end

    options = options_1;
  
    
