function [options] = process_options_structure(default_options, options, required_options, process_substructures)
% Return a structure the same as the structure options but with fields added from the structure default_options if they are not in options. Produces an error if any string in the cell array of strings required_options is not a field name in options. The second and third arguments are optional. Applies process_options_structure to fields of default_options that are structures if process_substructures = true.
% 
% Dependencies:
% 
% 2012-10-31 tebuck: Created.
% 2013-03-10 tebuck: Added optional third argument required_options.
% 2013-03-11 tebuck: Added optional fourth argument process_substructures.

  if nargin == 1
    options = default_options; 
    return
  end
  if ~exist('required_options', 'var')
    required_options = {}; 
  end
  if ~exist('process_substructures', 'var')
    process_substructures = false; 
  end
  
  option_names = fieldnames(default_options); 
  for index = 1:length(option_names)
    current_option = option_names{index}; 
    if ~isfield(options, current_option)
      if any(strcmp(current_option, required_options))
        error('Option %s is required', current_option)
      else
        options = setfield(options, current_option, getfield(default_options, current_option)); 
      end
    else
      % fprintf('isfield!\n')
      % keyboard
      if process_substructures && isfield(default_options, current_option) && isstruct(default_options.(current_option))
        % a = options.(current_option)
        % b = default_options.(current_option)
        options.(current_option) = process_options_structure(default_options.(current_option), options.(current_option), {}, process_substructures); 
        % b = default_options.(current_option)
        % a = options.(current_option)
        % keyboard
      end
    end
  end

    
