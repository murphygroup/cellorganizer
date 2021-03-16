function sbmlStruct = createImg2SBMLstruct(images,param)
%CREATEIMG2SBMLSTRUCT This function returns a structure for saving as SBML-Spatial
%
%Inputs:
%images = cell array containing images to be saved as SBML-spatial
%
%Outputs:
%framework = a struct for framework objects with triangulated mesh type
%encoding
%instance2SBML.m to generate a SBML-Spatial instance.
%

%Author: Devin Sullivan April 25, 2013
%
% Copyright (C) 2014 Murphy Lab
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



if nargin<2
    warning('No compartment names provided! Naming with sequential numbers')
end
    
%Check ordinality
if ~isfield(param,'ordinal')
    warning('No ordinality provided. Assigning all ordinals to 0')
    param.ordinal = [];
end


%Set up framework
sbmlStruct = struct;
sbmlStruct.name = 'frameworkMesh';

sbmlStruct.list = struct;
for i = 1:length(images)
    sbmlStruct.list(i).type = 'triangle mesh';
    sbmlStruct.list(i).name = param.SBML_Name{i};
    sbmlStruct.list(i).img = images{i};
     sbmlStruct.list(i).class = param.SBML_Name{i};
    if ~isfield(param,ordinal) || isempty(param.ordinal)
        sbmlStruct.list(i).ordinal = 0;
    else
        sbmlStruct.list(i).ordinal = param.ordinal{i};
    end
end

end