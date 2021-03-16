function answer = check_models_dimensionality( files )
% CHECK_MODELS_DIMENSIONALITY Helper function that returns true if all
% models in the files cell array have the same dimensionality. Otherwise,
% false.

% icaoberg@cmu.edu
%
% Copyright (C) 2018 Murphy Lab
% Computational Biology Department
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

answer = false;
if ~isa( files, 'cell' )
    warning( 'Input argument files must be a cell array' );
    return
end

if isempty( files )
    warning( 'Input argument files cannot be empty' );
    return
end

try
    load( files{1} );
    dimensionality = model.dimensionality;
    clear model
catch
    warning( 'Unable to extract dimensionality' );
    return
end

for i=2:1:length(files)
    load( files{i} )
    if ~strcmpi( dimensionality, model.dimensionality )
    end
end

answer = true;
return

end%check_models_dimensionality