function [ result ] = meshContainsPoint( mesh, position )
%MESHCONTAINSPOINT Test if mesh contains position.
%
% Inputs
% ------
% mesh     = a struct with field vertices and faces representing a closed, two-manifold mesh
% position = coordinates to test
%
% Outputs
% -------
% result = boolean indicating whether position is inside mesh


% Authors: Taraz Buck
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


[~, intersections_positions, ~] = intersectLineMesh3d([position, 1, 0, 0], mesh.vertices, mesh.faces);
% Odd number of intersections along a ray implies origin is within mesh
result = mod(sum(intersections_positions > 0), 2) == 1;


end

