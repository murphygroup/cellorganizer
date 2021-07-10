function model = param2model_3D( cell_params, options )
% PARAM2MODEL_3D Initial draft for param2model_3D.

% Eventually we would like an aribitrary model specified in options, and
% this function navigates the model structure and constructs it in the
% correct order.

% Copyright (C) 2016-2020 Murphy Lab
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

% May 25, 2016 I. Cao-Berg Included reference to model class and type
% Sep 17, 2018 Xiongtao Ruan Add 3D SPHARM-RPDM model

% options.if_skip_cell_nuclear_model=1;
if (isfield( options, 'nucleus') && ...
        isfield( options.nucleus, 'type') && ...
        strcmpi(options.nucleus.type, 'diffeomorphic')) || ...
        (isfield( options, 'cell') && ...
        isfield( options.cell, 'type') && ...
        strcmpi(options.cell.type, 'diffeomorphic'))
    isdiffeomorphic = true;
else
    isdiffeomorphic = false;
end

if (isfield( options, 'nucleus') && isfield( options.nucleus, 'type') && strcmpi(options.nucleus.type, 'spharm_rpdm')) ...
        || (isfield( options, 'cell')   && isfield( options.cell, 'type')    && strcmpi(options.cell.type, 'spharm_rpdm'))
    isrpdm = true;
else
    isrpdm = false;
end

if isempty(cell_params)
    warning(['Cell parameterization is empty. ' ...
        'Either images not found or unable to extract ' ...
        'parameters from images.'])
    model = [];
    return
end

if iscell(cell_params(1)) && ischar(cell_params{1})
    try
        cell_param = load(cell_params{1});
    catch
        warning('Unable to load parameters. Check the parameterization exist or is not empty.');
        cell_param = [];
        model = struct([]);
        return
    end
elseif iscell(cell_params(1))
    cell_param = cell_params{1};
else
    cell_param = cell_params(1);
end

%%Cell and Nuclear Shape Model
if isdiffeomorphic
    components = {};
    c = 1;
    if (isfield( options, 'cell') && ...
            isfield( options.cell, 'type') && ...
            strcmpi(options.cell.type, 'diffeomorphic'))
        components{c} = 'cell';
        c = c+1;
    end

    if (isfield( options, 'nucleus') && ...
            isfield( options.nucleus, 'type') && ...
            strcmpi(options.nucleus.type, 'diffeomorphic'))
        components{c} = 'nuc';
        c = c+1;
    end

    options_img = struct('components', {components}, 'format', 'indexed');
    options.imfunc = @(x) param2img(cell_params{x}, options_img);

    model.cellShapeModel = train_diffeomorphic_model(options);
    model.nuclearShapeModel = model.cellShapeModel;
elseif isrpdm
    components = {};
    if (isfield( options, 'cell')   && isfield( options.cell, 'type')    && strcmpi(options.cell.type, 'spharm_rpdm'))
        components{end + 1} = 'cell';
    end
    if (isfield( options, 'nucleus') && isfield( options.nucleus, 'type') && strcmpi(options.nucleus.type, 'spharm_rpdm'))
        components{end + 1} = 'nuc';
    end
    options.components = components;

    model.cellShapeModel = train_spharm_rpdm_model(cell_params, options);

    if ismember( options.train.flag, {'nuclear', 'framework', 'all'} )
    	model.nuclearShapeModel = model.cellShapeModel;
    	model.nuclearShapeModel.class = options.nucleus.class;
    	model.nuclearShapeModel.type = options.nucleus.type;
    	model.nuclearShapeModel.resolution = options.model.resolution;

	if isempty( model.nuclearShapeModel.name )
    		model.nuclearShapeModel.name = 'UNSET';
   	end

    	model.nuclearShapeModel.id = options.nucleus.id;
	if isempty( model.nuclearShapeModel.id )
    		model.nuclearShapeModel.id = uuidgen();
    	end
    end

    if ismember( options.train.flag, {'cell', 'framework', 'all'} )
    	model.cellShapeModel.class = options.cell.class;
    	model.cellShapeModel.type = options.cell.type;
    	model.cellShapeModel.resolution = options.model.resolution;

	if isempty( model.cellShapeModel.name )
        	model.cellShapeModel.name = 'UNSET';
        end

        model.cellShapeModel.id = options.cell.id;
        if isempty( model.cellShapeModel.id )
                model.cellShapeModel.id = uuidgen();
        end
    end
% elseif isfield(options, 'model') && isfield(options.model, 'name') && strcmp(options.model.name,'spharm_obj')
%     disp('generating spharm obj model');
%     model.spharm_obj_model = spharm_obj_model(options);
elseif isfield( options, 'if_skip_cell_nuclear_model' ) && options.if_skip_cell_nuclear_model
    disp('Skip creation of cell/nuclear model');
    model = {};
else
    if isfield(options, 'nucleus')
        nuc_model_save = [options.paramdir filesep 'nuc.mat'];
        if ~exist(nuc_model_save, 'file')
            switch options.nucleus.type
                case 'cylindrical_surface'
                    nuclear_model = param2cylsurf(cell_params, options);
                otherwise
                    warning(['Unsupported cell model type ' options.nucleus.type '. Returning empty model.'])
                    nuclear_model = [];
            end
            save(nuc_model_save, 'nuclear_model')
        else
            load(nuc_model_save)
        end

        model.nuclearShapeModel = nuclear_model;
        model.nuclearShapeModel.class = options.nucleus.class;
        model.nuclearShapeModel.type = options.nucleus.type;
        model.nuclearShapeModel.resolution = cell_param.options.model.resolution;
    end

    if isfield(options, 'cell')
        cell_model_save = [options.paramdir filesep 'cell.mat'];
        if ~exist(cell_model_save, 'file')
            switch options.cell.type
                case 'ratio'
                    cell_model = param2cell_ratio_model(cell_params, options);
                otherwise
                    warning(['Unsupported cell model type ' options.cell.type '. Returning empty model.'])
                    cell_model = [];
            end

            save(cell_model_save, 'cell_model')
        else
            load(cell_model_save)
        end
        model.cellShapeModel = cell_model;
        model.cellShapeModel.class = options.cell.class;
        model.cellShapeModel.type = options.cell.type;
        model.cellShapeModel.resolution = cell_param.options.model.resolution;
    end
end

%%Protein model
%icaoberg 01/02/2016
%first check that we are going to be using the protein pattern
if isfield(options, 'train') && isfield(options.train, 'flag') && ...
        (strcmpi( options.train.flag, 'all' )||strcmpi( options.train.flag, 'protein' ))
    %second verify the field is not empty
    if isfield(options, 'protein') && isfield(options.protein, 'type') && ...
            ~isempty(options.protein.type)
        protein_model_save = [options.paramdir filesep options.protein.type '.mat'];
        if ~exist(protein_model_save, 'file')
            switch options.protein.type
                case 'spharm_obj'
                    protein_model = spharm_obj_model(options);
                case 'gmm'
                    protein_model = param2protein_gaussian_model_3D(cell_params, options);
                case 'microtubule_growth'
                    protein_model = param2microtubule_model_3D(cell_params, options);
                otherwise
                    warning(['Unsupported protein model type ' options.protein.type '. Returning empty model.'])
                    protein_model = [];
            end

            model_type = options.protein.type;
            save(protein_model_save, 'protein_model', 'model_type')
        else
            load(protein_model_save)
        end

        model.proteinModel = protein_model;
        model.proteinModel.dimensionality = '3D';
        model.proteinModel.class = options.protein.class;
        model.proteinModel.type = options.protein.type;
        model.proteinModel.resolution = cell_param.options.model.resolution;
    end
end
end
