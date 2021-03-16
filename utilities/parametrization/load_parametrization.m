function parametrization = load_parametrization( parametrization )
%LOAD_PARAMETRIZATION Loads and returns the parametrization

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
end%load_parametrization
