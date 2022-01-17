function [ result ] = instance2SBML3_mod( CSGdata, meshData, models, imgs, savepath, SBMLfile, options )
%INSTANCE2SBML3_MOD Writes data struct to a SBML-Spatial file. If an existing file is provided
%this method will append to the end of the file.
%
% Inputs
% ------
% CSGdata = struct created using createSBMLstruct.m in slml2img.m
% Meshdata = cell array containing image arrays for mesh type objects to be saved
% savepath = output filename
% SBMLfile = optional path to existing file to append to.
%
% Outputs
% -------
% result = boolean flag indicating success
% saves .xml file in SBML-Spatial format.

% Authors: Rohan Arepally and Devin Sullivan
%
% Copyright (C) 2012-2017 Murphy Lab
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

%9/18/13 D. Sullivan - formatting and support for parametric objects
%9/24/13 D. Sullivan - edited for VCell compatability.

%D. Sullivan 9/10/14 - This is a flag to determine if we will create
%compartments in the xml file. If false, no compartments are created - this
%is not recommended for cases when using the model for simulation. If there
%is an SBML file already being appended to it will have the appropriate
%compartments.


debug = false;
% debug = true;

if ~debug
    warning('off', 'CellOrganizer:instance2SBML3_mod');
end


addCompartments = 0;

%D. Sullivan 9/10/14 - how many dimensions do you have?
ndim = 3;

if nargin==0
    warning('No input arguments given to instance2SBML3_mod. Nothing to do.');
    return
elseif nargin==1
    warning('No savepath given to instance2SBML3_mod. Defaulting to "./model.xml"');
    savepath = './model.xml';
    SBMLfile = 1;
end

if nargin<4
    SBMLfile = 1;
end

if isempty(SBMLfile)
    SBMLfile = 1;
end

%Alias s to 'spatial:'. this will be used to denote spatial specific
%attributes and nodes.
s = 'spatial:';
% listedCompartments = 0;
param.prefix = s;
param.ndim = 3;

if ~isempty(CSGdata) && ~isempty(meshData)
    param.mixed = 1;
end


% Names of loaded models other than cell and nucleus
model_names = cell(length(models), 1);
model_cytonuclearflags = cell(length(models), 1);
for i = 1:length(models)
    model = models{i};
    % name = ['model', num2str(i), '_', strrep(model.filename, '.', '_')];
    if isfield(model, 'filename')
        % Assume no duplicated filenames, make insensitive to input order
        name = strrep(model.filename, '.', '_');
    elseif isfield(model.documentation, 'original_files')
        name = [];
        for j = 1:length(model.documentation.original_files)
            [original_file_path, original_file_name, original_file_ext] = fileparts(model.documentation.original_files{j});
            if j > 1
                name = [name, '_'];
            end
            name = [name, original_file_name, '_', strrep(original_file_ext, '.', '')];
        end
    else
        name = ['model', num2str(i-1)];
    end
    cytonuclearflag = 'all';
    if isfield(model, 'proteinModel') && isfield(model.proteinModel, 'cytonuclearflag')
        cytonuclearflag = model.proteinModel.cytonuclearflag;
    end
    model_names{i} = name;
    model_cytonuclearflags{i} = cytonuclearflag;
end
options.output.SBMLModelNames = model_names;
options.output.SBMLModelCytonuclearflags = model_cytonuclearflags;


if ~isempty(CSGdata)
    options.output.use_individual_meshes = false;
    % options.output.use_individual_meshes = true;
    options.output.model_names = model_names;
    [CSGdata, meshData] = convertCSGToMesh(CSGdata, meshData, models, imgs, options);
    CSGdata = struct();
end


% Translate CSGdata and meshData names if user provide options.output.SBMLTranslations
translations = options.output.SBMLTranslations;
name_map = containers.Map('KeyType', 'char', 'ValueType', 'char');
% Defaults
name_map('EC') = 'EC';
name_map('cell') = 'cell';
name_map('nucleus') = 'nucleus';
for i = 1:length(models)
    name_map(model_names{i}) = model_names{i};
end
if ~isempty(translations)
    % Translations in options
    for i = 1:size(translations, 1)
        name_map(translations{i, 1}) = translations{i, 2};
    end
    if ~isempty(CSGdata)
        CSGdata_fieldnames = fieldnames(CSGdata);
        for i = 1:length(CSGdata_fieldnames)
            CSGdata_fieldname = CSGdata_fieldnames{i};
            CSGdata_item = CSGdata.(CSGdata_fieldname);
            CSGdata_list = CSGdata_item.list;
            for j = 1:length(CSGdata_list)
                CSGdata_list_item = CSGdata_list(j);
                CSGdata_list_item.name = translateWithDefaultIdentity(name_map, CSGdata_list_item.name);
                CSGdata_list(j) = CSGdata_list_item;
            end
            CSGdata_item.list = CSGdata_list;
            CSGdata = rmfield(CSGdata, CSGdata_fieldname);
            CSGdata_fieldname = translateWithDefaultIdentity(name_map, CSGdata_fieldname);
            CSGdata.(CSGdata_fieldname) = CSGdata_item;
        end
    end
    if ~isempty(meshData)
        for i = 1:length(meshData)
            single_meshData = meshData(i);
            single_meshData.name = translateWithDefaultIdentity(name_map, single_meshData.name);
            single_meshData_list = single_meshData.list;
            for j = 1:length(single_meshData_list)
                meshData_list_item = single_meshData_list(j);
                meshData_list_item.name = translateWithDefaultIdentity(name_map, meshData_list_item.name);
                single_meshData_list(j) = meshData_list_item;
            end
            single_meshData.list = single_meshData_list;
            meshData(i) = single_meshData;
        end
    end
end
options.output.SBMLNameMap = name_map;


[docNode,docRootNode,wrapperNode,GeowrapperNode,options] = setupSBML3(CSGdata,meshData,models,imgs, SBMLfile,param,options);


%%%Define the Geometry
geometryDefNode = docNode.createElement([param.prefix, 'listOfGeometryDefinitions']);
GeowrapperNode.appendChild(geometryDefNode);
% Treat all cases as mixed
geometryDefNodeMixed = docNode.createElement([param.prefix, 'mixedGeometry']);
geometryDefNodeMixed.setAttribute([s,'id'],['id' sprintf('%03d',num2str(randi([0 10000000]),1))]);
% Does not pass SBML validator
if options.output.SBMLSpatialVCellCompatible
    geometryDefNodeMixed.setAttribute('id',geometryDefNodeMixed.getAttribute([s,'id']))
end;
geometryDefNodeMixed.setAttribute([s,'isActive'],'true');
geometryDefNode.appendChild(geometryDefNodeMixed);
geometryDefNode = docNode.createElement([param.prefix, 'listOfGeometryDefinitions']);
geometryDefNodeMixed.appendChild(geometryDefNode);

ordinals = [];
domainlist = {};
objs = {};

%%%Set up the ListOfDomainTypes node
ListOfDomainTypesNode = docNode.createElement([s,'listOfDomainTypes']);
GeowrapperNode.appendChild(ListOfDomainTypesNode);
ListOfDomains = docNode.createElement([s,'listOfDomains']);
GeowrapperNode.appendChild(ListOfDomains);

if ~isempty(CSGdata) && ~options.output.SBMLSpatialImage
    [docNode,GeowrapperNode,geometryDefNode,wrapperNode] = addCSGObjects3(CSGdata,docNode,geometryDefNode,GeowrapperNode,wrapperNode,options);
    %Add surface area and volume for the class of the same "name" aka domainType
    
    %grab the ordinals
    domainlist = fieldnames(CSGdata);
    domainlist2 = {};
    for i = 1:length(domainlist)
        if ~strcmp(domainlist{i}, 'primitiveOnly')
            domainlist2{end+1} = domainlist{i};
        end
    end
    domainlist = domainlist2;
    for i = 1:length(domainlist)
        ordinals(i) = CSGdata.(domainlist{i}).ordinal;
        objs{i} = CSGdata.(domainlist{i}).list;
        %[docNode,wrapperNode] = addVol_SurfArea(docNode,wrapperNode,...
        %domainlist{i},CSGdata.(domainlist{i}).totvol,CSGdata.(domainlist{i}).totsa);
    end
end

%D. Sullivan 9/10/14  - added "primitiveOnly" for fast computations here
if (~isfield(CSGdata,'primitiveOnly') || CSGdata.primitiveOnly==0) && ...
        ~isempty(meshData)
    if options.output.SBMLSpatialImage
        % Write image data (for import into Virtual Cell):
        [docNode,GeowrapperNode,geometryDefNode,wrapperNode] = ...
            addImageObjects3(CSGdata,meshData,docNode,models,imgs,geometryDefNode, ...
            GeowrapperNode,wrapperNode, options);
    else
        % Write mesh data (doesn't work when imported into Virtual Cell):
        [docNode,GeowrapperNode,geometryDefNode,wrapperNode] = ...
            addMeshObjects3(meshData,docNode,geometryDefNode, ...
            GeowrapperNode,wrapperNode, options);
    end
    
    warning('CellOrganizer:instance2SBML3_mod','Assumes meshData(1) is frameworkMesh');
    for i = 1:length(meshData(1).list)
        domainlist{end+1} = meshData(1).list(i).name;
        ordinals(end+1) = meshData(1).list(i).ordinal;
        objs{end+1} = meshData(1).list(i);
    end
end
%%%

%%%List of adjacent domains
%     ordinals = extractfield(CSGdata.(domainlist).list,'ordinal');
[ordinal_list,order] = sort(unique(ordinals));
ListOfAdjacentDomains = docNode.createElement([s,'listOfAdjacentDomains']);
surface_list = options.output.SBMLSurfaceList;
surfaces_compartment_pairs = options.output.SBMLSurfacesCompartmentPairs;
volumes_surfaces = options.output.SBMLVolumesSurfaces;
%{
for j=1:length(ordinal_list)-1
    startobjs = {objs{find(ordinals==ordinal_list(j))}};
    linkobjs = {objs{find(ordinals==ordinal_list(j+1))}};
    %ugh, this is going to be really slow, but I can't think of a
    %better way of coding it right now.
    for k = 1:length(startobjs)
        %This is for when we have more than one object in a domain for
        %example. It is very slow and messy, but I don't have time to
        %restructure everything right now.
        for kk = 1:length(startobjs{k})
            %         start_object = CSGdata.list(startobjs(k));
            start_name = [startobjs{k}(kk).name,num2str(kk)];
            %         start_type = startobjs(1).type;
            
            for i = 1:length(linkobjs)
                %See above note for why this is necessary. Sorry it's so slow
                %and shitty.
                for ii = 1:length(linkobjs{i})
                    %             link_object = CSGdata.list(linkobjs(i));
                    link_name = [linkobjs{i}(ii).name,num2str(ii)];
                    %             link_type = link_object(i).type;
                    AdjacentDomainNode = docNode.createElement([s,'adjacentDomains']);
                    AdjacentDomainNode.setAttribute([s,'id'],[start_name,'_',link_name]);
                    % Does not pass SBML validator
                    if options.output.SBMLSpatialVCellCompatible
                        AdjacentDomainNode.setAttribute('id',AdjacentDomainNode.getAttribute([s,'id']))
                    end;
                    AdjacentDomainNode.setAttribute([s,'domain1'],start_name);
                    AdjacentDomainNode.setAttribute([s,'domain2'],link_name);
                    ListOfAdjacentDomains.appendChild(AdjacentDomainNode);
                end
            end
        end
    end
    
end
%}

for i = 1:length(surface_list)
    surface_name = surface_list{i};
    surface_compartment_pair = surfaces_compartment_pairs(i, :);
    model_name_mapped = translateWithDefaultIdentity(name_map, surface_name);
    compartment_outside_name = surface_compartment_pair{1};
    compartment_inside_name = surface_compartment_pair{2};
    compartment_outside_name = translateWithDefaultIdentity(name_map, compartment_outside_name);
    
    % Add to listOfDomainTypes:
    DomainTypeNode = docNode.createElement([s,'domainType']);
    % DomainTypeNodeID = [name,'_domainType'];
    DomainTypeNodeID = [surface_name,'_domainType'];
    DomainTypeNode.setAttribute([s,'id'],DomainTypeNodeID);
    % Does not pass SBML validator
    if options.output.SBMLSpatialVCellCompatible
        DomainTypeNode.setAttribute('id',DomainTypeNode.getAttribute([s,'id']))
    end;
    DomainTypeNode.setAttribute([s,'spatialDimensions'],'2');
    ListOfDomainTypesNode.appendChild(DomainTypeNode);

    % Add to listOfDomains:
    DomainNode = docNode.createElement([s,'domain']);
    % DomainNode.setAttribute([s,'id'],[name,'_domain']);
    DomainNode.setAttribute([s,'id'],surface_name);
    % Does not pass SBML validator
    if options.output.SBMLSpatialVCellCompatible
        DomainNode.setAttribute('id',DomainNode.getAttribute([s,'id']))
    end;
    DomainNode.setAttribute([s,'domainType'],DomainTypeNodeID);
    % ListOfInteriorPoints = docNode.createElement([s,'listOfInteriorPoints']);
    % InteriorPoint = docNode.createElement([s,'interiorPoint']);
    ListOfDomains.appendChild(DomainNode);
    
    
    %{
    ListOfInteriorPointsNode = docNode.createElement([s,'listOfInteriorPoints']);
    DomainNode.appendChild(ListOfInteriorPointsNode);
    InteriorPointNode = docNode.createElement([s,'interiorPoint']);
    ListOfInteriorPointsNode.appendChild(InteriorPointNode);
    [ipi, ipj, ipk] = find(object.img, 1, 'first');
    InteriorPointNode.setAttribute([s,'coord1'],num2str(ipj * resolution_ijk(2)));
    InteriorPointNode.setAttribute([s,'coord2'],num2str(ipi * resolution_ijk(1)));
    InteriorPointNode.setAttribute([s,'coord3'],num2str(ipk * resolution_ijk(3)));
    %}
end

for i = 1:size(surfaces_compartment_pairs, 1)
    surface_compartment_pair = surfaces_compartment_pairs(i, :);
    surface_name = surface_list{i};
    
    AdjacentDomainNode = docNode.createElement([s,'adjacentDomains']);
    AdjacentDomainNode.setAttribute([s,'id'],[surface_compartment_pair{1},'__',surface_name]);
    AdjacentDomainNode.setAttribute([s,'domain1'],surface_compartment_pair{1});
    AdjacentDomainNode.setAttribute([s,'domain2'],surface_name);
    ListOfAdjacentDomains.appendChild(AdjacentDomainNode);
    
    AdjacentDomainNode = docNode.createElement([s,'adjacentDomains']);
    AdjacentDomainNode.setAttribute([s,'id'],[surface_name,'__',surface_compartment_pair{2}]);
    AdjacentDomainNode.setAttribute([s,'domain1'],surface_name);
    AdjacentDomainNode.setAttribute([s,'domain2'],surface_compartment_pair{2});
    ListOfAdjacentDomains.appendChild(AdjacentDomainNode);
end

if ListOfAdjacentDomains.getChildNodes().getLength() > 0
    GeowrapperNode.appendChild(ListOfAdjacentDomains);
end

% Simulation parameters
wrapper_annotation = docNode.createElement(['annotation']);
cb = 'cellblender:';
wrapper_annotation_gsc = docNode.createElement([cb,'simulationParameters']);
docRootNode.setAttribute(['xmlns:','cellblender'],'http://mcell.org/CellBlender');
wrapper_annotation_gsc.setAttribute([cb,'endTime'],num2str(options.output.SBMLEndTime));
wrapper_annotation_gsc.setAttribute([cb,'timeStep'],num2str(options.output.SBMLDefaultTimeStep));
wrapper_annotation_gsc.setAttribute([cb,'maxTimeStep'],num2str(options.output.SBMLMaxTimeStep));
wrapper_annotation_gsc.setAttribute([cb,'outputTimeStep'],num2str(options.output.SBMLOutputTimeStep));
wrapper_annotation.appendChild(wrapper_annotation_gsc);
wrapperNode.appendChild(wrapper_annotation);



%%%
%%%Append all the children nodes and save the file

wrapperNode.appendChild(GeowrapperNode);

% Rearrange spatial:geometry children
elements_to_rearrange = {};
elements_to_rearrange{end+1} = 'listOfAdjacentDomains';
elements_to_rearrange{end+1} = 'listOfCoordinateComponents';
elements_to_rearrange{end+1} = 'listOfDomains';
elements_to_rearrange{end+1} = 'listOfDomainTypes';
elements_to_rearrange{end+1} = 'listOfGeometryDefinitions';
elements_to_rearrange{end+1} = 'listOfSampledFields';
for i = 1:length(elements_to_rearrange)
    element_name = [s, elements_to_rearrange{i}];
    element = GeowrapperNode.getElementsByTagName(element_name);
    if element.getLength() > 0
        element = element.item(0);
        GeowrapperNode.appendChild(GeowrapperNode.removeChild(element));
    end
end

% Rearrange model children
elements_to_rearrange = {};
elements_to_rearrange{end+1} = 'listOfUnitDefinitions';
elements_to_rearrange{end+1} = 'listOfCompartments';
elements_to_rearrange{end+1} = 'listOfParameters';
for i = 1:length(elements_to_rearrange)
    element_name = elements_to_rearrange{i};
    element = wrapperNode.getElementsByTagName(element_name);
    if element.getLength() > 0
        element = element.item(0);
        wrapperNode.appendChild(wrapperNode.removeChild(element));
    end
end

% docRootNode.appendChild(GeowrapperNode);
docRootNode.appendChild(wrapperNode);
xmlwrite(savepath, docNode);
zip([savepath, '.zip'], savepath);

result = 1;
end

