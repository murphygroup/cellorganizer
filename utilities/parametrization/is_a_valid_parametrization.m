function answer = is_a_valid_parametrization( parametrization )
%IS_A_VALID_PARAMETRIZATION Returns true if it is a valid parametrization structure. False
%otherwise.

% Author: Ivan E. Cao-Berg (icaoberg@cmu.edu)
% Created: October 12, 2015
%
% Copyright (C) 2015 Murphy Lab
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
%
% For additional information visit http://murphylab.web.cmu.edu or
% send email to murphy@cmu.edu

% List of steps
% * step0: check input arguments
% * step1: check paramtare valid
% * step3: check name is a member of list

%% step0: check input arguments
answer = false;
if nargin ~= 1
    warning( 'Wrong number of input arguments.' );
    return
end

if isa( parametrization, 'char' )
    %check if this string points to a file on disk
    if exist( parametrization, 'file' )
        %if file exists, then try load it
        try
            load( parametrization )
            parametrization = cell_param;
            clear cell_param
        catch e
            warning( ['Unable to load: ' parametrization ] );
            getReport( e )
            return
        end
    else
        %when the string refers to a file on disk but the file does not
        %exist
        warning( ['Unable to find: ' parametrization ] );
        return
    end
elseif isstruct( parametrization )
    if isempty( parametrization )
        warning('Input argument parametrization is a structure but it is empty.');
        return
    end
else
    warning( 'Input argument parametrization must be either a string or a structure' );
    return
end

%check dimensionality
answer = check_dimensionality( parametrization );
if ~answer
    warning( 'Invalid parametrization due to dimensionality.' );
    return
end

%check nuclear component
answer = check_nuclear_component( parametrization );
if ~answer
    warning( 'Invalid parametrization due to nuclear component.' );
    return
end

%check cell component
answer = check_cell_component( parametrization );
if ~answer
    warning( 'Invalid parametrization due to cell model component.' );
    return
end

%check protein component
answer = check_protein_component( parametrization );
if ~answer
    warning( 'Invalid parametrization due to protein component.' );
    return
end

end%is_a_valid_parametrization

function answer = check_dimensionality( parametrization )
%Helper function that checks the dimensionality in the parametrization
answer = false;

if ~isfield( parametrization, 'dimensionality' )
    warning( 'Parametrization does not contain dimensationality' );
    return
end

dimensionalities = { '2D', '3D', '4D' };
answer = ismember( parametrization.dimensionality, ...
    dimensionalities );
if ~answer
    warning( 'Parametrization contained an invalid value for dimensionality' );
end
return

end%check_dimensionality

function answer = check_nuclear_component( parametrization )
%Helper function that checks the nuclear component of the parametrization
answer = false;

%check nuclear field exists
if ~isfield( parametrization, 'nucleus' )
    warning( 'Parametrization does not contain nuclear component' );
    return
end

%check nuclear model class field exists and we have a valid value
if ~isfield( parametrization.nucleus, 'class' )
    warning( 'Nuclear component does not contain model class' );
    return
else
    answer = is_a_model_class( parametrization.nucleus.class );
    if ~answer
        warning( 'Invalid nuclear model class value' );
        return
    end
end

%check nuclear model type field exists and we have a valid value
if ~isfield( parametrization.nucleus, 'type' )
    warning( 'Nuclear component does not contain model type' );
    return
else
    answer = is_a_model_type( parametrization.nucleus.type );
    if ~answer
        warning( 'Invalid nuclear model type value' );
        return
    end
end
end%check_nuclear_component

function answer = check_cell_component( parametrization )
%Helper function that checks the cell component of the parametrization
answer = false;

%check nuclear field exists
if ~isfield( parametrization, 'cell' )
    warning( 'Parametrization does not contain cell component' );
    return
end

%check cell model class field exists and we have a valid value
if ~isfield( parametrization.cell, 'class' )
    warning( 'Cellular component does not contain model class' );
    return
else
    answer = is_a_model_class( parametrization.cell.class );
    if ~answer
        warning( 'Invalid cell model class value' );
        return
    end
end

%check cell model type field exists and we have a valid value
if ~isfield( parametrization.cell, 'type' )
    warning( 'Cellular component does not contain model type' );
    return
else
    answer = is_a_model_type( parametrization.cell.type );
    if ~answer
        warning( 'Invalid cell model type value' );
        return
    end
end
end%check_cell_component


function answer = check_protein_component( parametrization )
%Helper function that checks the protein component of the parametrization
answer = false;

%check protein field exists
if ~isfield( parametrization, 'protein' )
    warning( 'Parametrization does not contain protein component' );
    return
end

%check protein model class field exists and we have a valid value
if ~isfield( parametrization.protein, 'class' )
    warning( 'Protein component does not contain model class' );
    return
else
    answer = is_a_model_class( parametrization.protein.class );
    if ~answer
        warning( 'Invalid protein model class value' );
        return
    end
end

%check protein model type field exists and we have a valid value
if ~isfield( parametrization.cell, 'type' )
    warning( 'Protein component does not contain model type' );
    return
else
    answer = is_a_model_type( parametrization.protein.type );
    if ~answer
        warning( 'Invalid protein model type value' );
        return
    end
end
end%check_protein_component
