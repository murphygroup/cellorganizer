function answer = is_a_model_type( name )
%IS_A_MODEL_TYPE Returns true if name is a valid model type name. False
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
% * step1: check options are valid
% * step3: check name is a member of list

%% step0: check input arguments
answer = false;
if nargin ~= 1
    warning( 'Wrong number of input arguments.' );
    return
end

%% step1: check options are valid
if ~isa( name, 'char' )
    warning( 'Input argument ''name'' must be a string.' );
    return
end

%% step3: check name is a member of list

%COMMENT_FOR_DEVELOPERS: when adding new types, add them to the list below

%list of valid model types
list_of_types = { 'gmm', 'medial_axis', 'ratio', 'cylindrical_surface', ...
    'diffeomorphic', 'microtubule_growth', ...
    'spatial_point_process_microtubule' };

answer = ismember( name, list_of_types );

end%is_a_model_type