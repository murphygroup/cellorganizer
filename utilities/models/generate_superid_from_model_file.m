function id = generate_superid_from_model_file(filename)
% GENERATE_SUPERID_FROM_MODEL_FILE Takes a filename as input argument, and returns a super string. 
%
% The super string should be composed from the
% model.id
% model.nuclearShapeModel.id
% model.cellShapeModel.id
% model.proteinShape.id
%
% List Of Input Arguments     Descriptions
% -----------------------     ------------
% filename                    Name of the file that will be used
%
% Author: Juan Pablo Hinojosa(jhinojos@andrew.cmu.edu), Rita Chen(ruijiac@andrew.cmu.edu)
%
% Copyright (C) 2007-2018 Murphy Lab
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

% load the model
S = load(filename);

% check whether cellShapeModel, nuclearShapeModel, proteinModel exist

if (isfield(S.model,'cellShapeModel') && isfield(S.model.cellShapeModel,'id') && (~isempty(S.model.cellShapeModel.id)))
    ID_cellShapeModel = S.model.cellShapeModel.id;
else
    ID_cellShapeModel = "";
end

if (isfield(S.model,'nuclearShapeModel') && isfield(S.model.nuclearShapeModel,'id') && (~isempty(S.model.nuclearShapeModel.id)))
    ID_nuclearShapeModel = S.model.cellShapeModel.id;
else
    ID_nuclearShapeModel = "";
end


if (isfield(S.model,'proteinModel') && isfield(S.model.proteinModel,'id') && (~isempty(S.model.proteinModel.id)))

    ID_proteinModel = S.model.cellShapeModel.id;
else
    ID_proteinModel = "";
end

id_list = ['m' S.model.id '-n' ID_nuclearShapeModel '-c' ID_cellShapeModel '-p' ID_proteinModel];

try
	id = strjoin(id_list,'');
catch
	id = id_list;
end
