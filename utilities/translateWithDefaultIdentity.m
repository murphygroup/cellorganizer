function value = translateWithDefaultIdentity(map_object, key)
%TRANSLATEWITHDEFAULTIDENTITY Convertes CellOrganizer geometry in CSGdata format to meshData format.
%
% Inputs
% ------
% map_object = 
% key        = 
%
% Outputs
% -------
% value = 


% Authors: Taraz Buck
%
% Copyright (C) 2019 Murphy Lab


if map_object.isKey(key)
    value = map_object(key);
else
    value = key;
end

end
