function model = param2model_2D( cell_params, options )
% PARAM2MODEL_2D Initial draft for param2model_2D.
%
% Eventually we would like an aribitrary model specified in options, and
% this function navigates the model structure and constructs it in the
% correct order.

% Copyright (C) 2015-2017 Murphy Lab
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
% For additional information visit http://murphylab.web.cmu.edu/ or
% send email to murphy@cmu.edu
%
% Feb 17, 2018 X. Ruan add code for 2D PCA model
% 8/14/2022 R.F.Murphy add check for empty file list

disp(upper('Checking if we are training a diffeomorphic model'));
if (fieldisequal( options, 'nucleus.type', 'diffeomorphic', false)) ...
        || (fieldisequal( options, 'cell.type', 'diffeomorphic', false))
    isdiffeomorphic = true;
    disp('Training diffeomorphic model');
else
    isdiffeomorphic = false;
    disp('Not training a diffeomorphic model');
end

disp(upper('Checking if we are training a PCA model'))
if (fieldisequal( options, 'nucleus.type', 'pca', false)) ...
        || (fieldisequal( options, 'cell.type', 'pca', false))
    ispca = true;
    disp('Training a PCA model')
else
    disp('Not training a PCA model');
    ispca = false;
end

model.dimensionality = options.dimensionality;

if isempty(cell_params)
    model = [];
    warning('Unable to extract parameters from input file list')
    return
end

if isempty(cell_params)
    model = [];
    warning('Unable to extract parameters')
    return
end

%loop through parameterizations until we find one that exists
c = 1;
breakflag = false;
while ~breakflag
    if iscell(cell_params(c)) && ischar(cell_params{c}) && exist(cell_params{c}, 'file')
        cell_param = load(cell_params{c}, 'options');
        breakflag = true;
    elseif iscell(cell_params(c)) && isstruct(cell_params{c})
        cell_param = cell_params{c};
        breakflag = true;
    elseif isstruct(cell_params(c))
        cell_param = cell_params(c);
        breakflag = true;
    end
    c = c+1;
end

disp('Building generative models from parameterization');
if isdiffeomorphic
    model = build_diffeomorphic_model( cell_params, cell_param, options );
elseif ispca
    model = build_pca_model( cell_params, cell_param, options );
else
    model = build_other_model( cell_params, cell_param, options );
end
end%param2model_2D

function model = build_diffeomorphic_model( cell_params, cell_param, options )
components = {};

if (fieldisequal( options, 'cell.type', 'diffeomorphic', false))
    % xruan 01/05/2016 change components{end} = 'cell'; to components{end + 1} = 'cell';
    components{end + 1} = 'cell';
end

if (fieldisequal( options, 'nucleus.type', 'diffeomorphic', false))
    % xruan 01/05/2016 change components{end} = 'nucleus'; to components{end + 1} = 'nucleus';
    components{end + 1} = 'nuc';
end

% xruan 01/05/2016 change component to components to unify the name
% make options_img as a single structure
options_img = struct('components', {components}, 'format', 'indexed');
options.imfunc = @(x) param2img(cell_params{x}, options_img);
% options.imfunc(1);

model.cellShapeModel = train_diffeomorphic_model(options);
model.nuclearShapeModel = model.cellShapeModel;
model.nuclearShapeModel.class = options.nucleus.class;
model.nuclearShapeModel.type = options.nucleus.type;
model.nuclearShapeModel.name = options.nucleus.name;
if isempty( model.nuclearShapeModel.name )
    model.nuclearShapeModel.name = 'UNSET';
end
model.nuclearShapeModel.id = options.nucleus.id;
if isempty( model.nuclearShapeModel.id )
    model.nuclearShapeModel.ID = uuidgen();
end
end%build_diffeomorphic_model

function model = build_pca_model( cell_params, cell_param, options )
components = {};
if (fieldisequal( options, 'cell.type', 'pca', false))
    components{end + 1} = 'cell';
end
if (fieldisequal( options, 'nucleus.type', 'pca', false))
    components{end + 1} = 'nuc';
end
options.components = components;

model.cellShapeModel = train_pca_2d_model(cell_params, options);
model.nuclearShapeModel = model.cellShapeModel;
model.nuclearShapeModel.class = options.nucleus.class;
model.nuclearShapeModel.type = options.nucleus.type;
model.cellShapeModel.class = options.cell.class;
model.cellShapeModel.type = options.cell.type;
if isempty( model.nuclearShapeModel.name )
    model.nuclearShapeModel.name = 'UNSET';
end
model.nuclearShapeModel.id = options.nucleus.id;
if isempty( model.nuclearShapeModel.id )
    model.nuclearShapeModel.id = uuidgen();
end
end%build_pca_model

function model = build_other_model( cell_params, cell_param, options )

if ismember( options.train.flag, {'nuclear', 'framework', 'all', 'protein'} )
    if isfield(options, 'nucleus')
        nuc_model_save = [options.paramdir filesep 'nuc.mat'];
        if ~exist(nuc_model_save, 'file')
            switch options.nucleus.type
                case 'medial_axis'
                    options.model = ml_initparam(options.model, struct('nucleus', []));
                    options.model.nucleus = ml_initparam(options.model.nucleus, ...
                        struct('modelname','mxp','constknots',{{0.5,0.5}}));
                    
                    param_tmp = loadFiles(cell_params, 'nuc');
                    param_tmp = param_tmp(~cellfun(@isempty, param_tmp));
                    
                    component_param = cellfun(@(x) x.nuc, param_tmp, 'UniformOutput', false);
                    
                    nuclear_model = ml_trainshapemodel2D( component_param, options.model.nucleus );
                    
                    %                     nuclear_model = param2cylsurf(cell_params, options);
                otherwise
                    warning(['Unsupported cell model type ' options.nucleus.type '. Returning empty model.'])
                    nuclear_model = [];
            end
            save(nuc_model_save, 'nuclear_model')
        else
            load(nuc_model_save)
        end
        
        model.nuclearShapeModel = nuclear_model;
        model.nuclearShapeModel.resolution = cell_param.options.model.resolution;
        model.nuclearShapeModel.class = options.nucleus.class;
        model.nuclearShapeModel.type = options.nucleus.type;
        if isempty( model.nuclearShapeModel.name )
            model.nuclearShapeModel.name = 'UNSET';
        end
        model.nuclearShapeModel.id = options.nucleus.id;
        if isempty( model.nuclearShapeModel.id )
            model.nuclearShapeModel.id = uuidgen();
        end
    end%if/else nucleus
end

if ismember( options.train.flag, {'cell', 'framework', 'all'} )
    if isfield(options, 'cell')
        cell_model_save = [options.paramdir filesep 'cell.mat'];
        if ~exist(cell_model_save, 'file')
            switch options.cell.type
                case 'ratio'
                    param_tmp = loadFiles(cell_params, 'cell');
                    param_tmp = param_tmp(~cellfun(@isempty, param_tmp));
                    
                    component_param = cellfun(@(x) x.cell, param_tmp, 'UniformOutput', false);
                    
                    cell_model = param2cell_ratio_model_2D(component_param, options);
                otherwise
                    warning(['Unsupported cell model type ' options.cell.type '. Returning empty model.'])
                    cell_model = [];
            end
            
            save(cell_model_save, 'cell_model')
        else
            load(cell_model_save)
        end
        
        model.cellShapeModel = cell_model;
        model.cellShapeModel.resolution = cell_param.options.model.resolution;
        model.cellShapeModel.class = options.cell.class;
        model.cellShapeModel.type = options.cell.type;
        model.cellShapeModel.type = options.nucleus.type;
        if isempty( model.cellShapeModel.name )
            model.cellShapeModel.name = 'UNSET';
        end
        model.cellShapeModel.id = options.nucleus.id;
        if isempty( model.cellShapeModel.id )
            model.cellShapeModel.id = uuidgen();
        end
    end%if/else cell
end

if ismember( options.train.flag, {'all', 'protein'} )
    if isfieldr(options, 'protein.type') && ...
            ~isempty( options.protein.type )
        protein_model_save = [options.paramdir filesep options.protein.type '.mat'];
        if ~exist(protein_model_save, 'file')
            
            switch options.protein.type
                case 'gmm'
                    param_tmp = loadFiles(cell_params, {'prot', 'cell', 'imsize'});
                    param_tmp = param_tmp(~cellfun(@isempty, param_tmp));
                    
                    component_param_prot = cellfun(@(x) x.prot, param_tmp, 'UniformOutput', false);
                    component_param_cell = cellfun(@(x) x.cell, param_tmp, 'UniformOutput', false);
                    imsizes = cellfun(@(x) x.imsize, param_tmp, 'UniformOutput', false);
                    imsizes = vertcat(imsizes{:});
                    
                    protein_model = param2protein_gaussian_model_2D(component_param_prot, component_param_cell, imsizes, options);
                case 'network'
                    warning(['Unsupported protein model type "network" is currently not supported. Returning empty model.'])
                    %                 protein_model = param2microtubule_model_3D(cell_params, options);
                otherwise
                    warning(['Unsupported protein model type ' options.protein.type '. Returning empty model.'])
                    protein_model = [];
            end
            save(protein_model_save, 'protein_model')
        else
            load(protein_model_save)
        end
        
        model.proteinModel = protein_model;
        model.proteinModel.dimensionality = '2D';
        model.proteinModel.resolution = cell_param.options.model.resolution;
        model.proteinModel.class = options.protein.class;
        model.proteinModel.type = options.protein.type;
        
        if isempty( model.nuclearShapeModel.name )
            model.proteinModel.name = 'UNSET';
        end
        model.proteinModel.id = options.nucleus.id;
        if isempty( model.proteinModel.id )
            model.proteinModel.id = uuidgen();
        end
    end
    
end
end%build_other_model
