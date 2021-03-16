function answer = is_classic_model( model )
% IS_CLASSIC_MODEL Helper function that returns true if model is a classic
% 3D model

% Ivan E. Cao-Berg
%
% Copyright (C) 2018 Murphy Lab
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

if strcmpi( model.dimensionality, '3D' )
    try
        if  strcmpi( model.nuclearShapeModel.class, 'nuclear_membrane' ) && ...
                strcmpi( model.nuclearShapeModel.type, 'cylindrical_surface' ) && ...
                strcmpi( model.cellShapeModel.class, 'cell_membrane' ) && ...
                strcmpi( model.cellShapeModel.type, 'ratio' ) && ...
                strcmpi( model.proteinModel.class, 'vesicle' ) && ...
                strcmpi( model.proteinModel.type, 'gmm' )
            answer = true;
        else
            answer = false;
        end
    catch
        answer = false;
    end
elseif strcmpi( model.dimensionality, '2D' )
    try
        if  strcmpi( model.nuclearShapeModel.class, 'nuclear_membrane' ) && ...
                strcmpi( model.nuclearShapeModel.type, 'medial_axis' ) && ...
                strcmpi( model.cellShapeModel.class, 'cell_membrane' ) && ...
                strcmpi( model.cellShapeModel.type, 'ratio' ) && ...
                strcmpi( model.proteinModel.class, 'vesicle' ) && ...
                strcmpi( model.proteinModel.type, 'gmm' )
            answer = true;
        else
            answer = false;
        end
    catch
        answer = false;
    end
else
    answer = false;
end
end%is_classic_model