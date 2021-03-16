function primitive = createSBMLstruct(model,currentmodel,modelclass,primitive)
%CREATESBMLSTRUCT This function returns a primitive structure for saving as SBML-Spatial
%
%Inputs: 
%model = Protein model struct from which the temp primitive were generated
%currentmodel = suffix which the primitive temp file contains
%distinguishing it from other models being run. 
%
%Outputs:
%primitive = a struct for primitive object types for input into
%instance2SBML.m to generate a SBML-Spatial instance.
%

%Author: Devin Sullivan September 16, 2013
%Edited: 
%D. Sullivan 9/18/13 - Cleaned up code and added support for multiple
%models
%
% Copyright (C) 2013 Murphy Lab
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




%first load temp results 
try 
    load([pwd filesep 'temp' filesep 'primitives' currentmodel '.mat']);
catch 
    error('Unable to locate temp primitive file for the specified model.')
end

if nargin<4    
    primitive = struct;
end

if nargin<3 || isempty(modelclass)
    modelclass = model.class
end

% if ~isfield(primitive,'class')
%     primitive.class = modelclass;
% end
% 
% if ~isfield(primitive,'name')
%     primitive.name = modelclass;%model.name;
% end
% if ~isfield(primitive,(modelclass))
%     primitive.(modelclass) = struct;
% end
% if ~isfield(primitive.(modelclass),'list')
%     primitive.(modelclass).list = struct;
% end
% if ~isfield(primitive.(modelclass).list,'type')
%     primitive.(modelclass).list.type = struct;
% end
% if ~isfield(primitive.(modelclass).list,'name')
%     primitive.(modelclass).list.name = struct;
% end
% if ~isfield(primitive.(modelclass).list,'position')
%     primitive.(modelclass).list.position = struct;
% end
% if ~isfield(primitive.(modelclass).list,'rotation')
%     primitive.(modelclass).list.rotation = struct;
% end
% if ~isfield(primitive.(modelclass).list,'scale')
%     primitive.(modelclass).list.scale = struct;
% end

% model.proteinModel.cytonuclearflag default is said to be 'all' in slml2img comments but 'cyto' in demo comments, and options.cytonuclearflag defaults to 'cyto' in get_cellorganizer_default_parameters.m
model = ml_initparam(model,struct('cytonuclearflag','cyto'));

%Determine the ordinal relative to which compartments objects are allowed
%ordinal is indexed at 0 for extracellular 
switch model.cytonuclearflag
    case 'cyto'
        %inside cell but not nuc
        ordinal = 2;
    case 'nuc'
        %inside nuc
        ordinal = 3;
    case 'all'
        %this doesn't really have a proper definition in SBML-Spatial
        %because it allows objects to directly overlap.
        ordinal = 1.5;
    otherwise
        warning('Unrecognized value for cytonuclearflag specified in protein model, setting ordinal to default of 2.');
        ordinal = 2;
end

%Loop through each object
% nprims = length(primitive.(modelclass).list.name);
currobj = 1;
cubicres = repmat(min(model.resolution),size(model.resolution,1),size(model.resolution,2));
objsizescaled = objsizevec.*repmat(cubicres,size(objsizevec,1),1);
objposscaled = objposvec.*repmat(cubicres,size(objsizevec,1),1);
for i = 1:size(objsizevec,1)
    
%     if ~isempty(primitive.(modelclass).list(i))
    %Recognize the primitive type to be used based on the model.type
    switch model.type
        case 'gmm'
            primitive.(modelclass).list(i).type = 'sphere';
        otherwise
            warning('non-recognized SBML model type. Using blank string.');
            primitive.(modelclass).list.type(i) = '';
    end

    primitive.(modelclass).list(i).name = [modelclass];
    
%     primitive.(modelclass).list(i).position = objposvec(currobj,:).*model.resolution;
    primitive.(modelclass).list(i).position = objposscaled(currobj,:);
    primitive.(modelclass).list(i).rotation = objrotvec(currobj,:);%This is scale invarient 
%     primitive.(modelclass).list(i).scale = objsizevec(currobj,:).*model.resolution;
    primitive.(modelclass).list(i).scale = objsizescaled(currobj,:);%[0.05,0.05,0.04];
    primitive.(modelclass).list(i).rotationmatrix = squeeze(objrotmat(currobj,:,:));
    primitive.(modelclass).list(i).covariancematrix = squeeze(objcovmat(currobj,:,:));
    
    
    %D. Sullivan - 11/10/14
    %Compute the volume and surface area 
    primitive.(modelclass).list(i).vol = 4/3*pi*prod(objsizescaled(currobj,:));
    %for SA we will use the Knud Thomsen's formula (error at most 1.061%)
    p = 1.6075;
    primitive.(modelclass).list(i).sa = 4*pi*(((objsizescaled(currobj,1)*objsizescaled(currobj,2))^p+...
        (objsizescaled(currobj,1)*objsizescaled(currobj,3))^p+...
        (objsizescaled(currobj,2)*objsizescaled(currobj,3))^p)/3)^(1/p);
    
    
    currobj = currobj+1;
end

 primitive.(modelclass).ordinal = ordinal;

primitive.(modelclass).totvol = sum(extractfield(primitive.(modelclass).list,'vol'));
primitive.(modelclass).totsa = sum(extractfield(primitive.(modelclass).list,'sa'));
