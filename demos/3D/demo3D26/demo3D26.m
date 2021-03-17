function answer = demo3D26()
% demo3D26
%
% This function displays a shape space of some dimensionality. This demo
% uses the model trained in Johnson 2015.
%
% Input 
% -----
% * a CellOrganizer diffeomorphic model
%
% Output
% ------
% * a display of the shape space

% Copyright (C) 2016-2017 Murphy Lab
% Computational Biology Department
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT MODIFY THIS BLOCK
if ~isdeployed()
  current_path = which(mfilename);
  [current_path, filename, extension] = fileparts( current_path );
  cd(current_path);
end

disp( 'demo3D26' );
disp( 'The estimated running time is 2 hours and 40 minutes. Please wait.' );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modelname = 'helamodel_9_16_15.mat';

disp(['Loading ' modelname '. This might take a few minutes']);
try
    load(['../../../models/3D/diffeomorphic/' modelname])
catch
    warning( 'Unable to load model. Exiting demo.' );
    getReport(err)
    return
end
disp([modelname ' loaded.'])

disp(['Warning: building the display of the shape space ' ...
    'requires considerable resources. Please be patient.']);
test_reconstruction_methods( model.cellShapeModel.distances_incomplete )
answer = true;
end%demo3D26
