function [docNode,docRootNode,wrapperNode,GeowrapperNode,options] = setupSBML3(CSGdata,meshData,models,imgs,SBMLfile,param,options)


debug = false;
% debug = true;

if ~debug
    warning('off', 'CellOrganizer:setupSBML3');
end


%D. Sullivan 9/10/14 - how many dimensions do you have?
%D. Sullivan 9/11/14 added param structure to make initialization easier
ndim = param.ndim;
s = param.prefix;

img_size = size(imgs{1});
warning('CellOrganizer:setupSBML3','Assumes meshData(1) is frameworkMesh');
resolution_ijk = options.resolution.cubic([2,1,3]);
img_size_real = img_size .* resolution_ijk;

model_names = options.output.SBMLModelNames;
model_cytonuclearflags = options.output.SBMLModelCytonuclearflags;
name_map = options.output.SBMLNameMap;

if(ischar(SBMLfile))
    docNode = xmlread(SBMLfile);
    docRootNode = docNode.getDocumentElement;
    docRootNode = XMLremoveWhitespaceNodes(docRootNode);
    allListitems = docNode.getElementsByTagName('model');
    wrapperNode = allListitems.item(0);
    allListitems = wrapperNode.getElementsByTagName('listOfCompartments');
    ListOfCompartments = allListitems.item(0);
    %If we've already defined the compartment list in the SBML file, no
    %need to re-do it.
    if length(ListOfCompartments~=0)
        listedCompartments = 1;
    else 
        listedCompartments = 0;
    end
    
    % Prevent Virtual Cell error
    if options.output.SBMLSpatialVCellCompatible
        speciesReferences = wrapperNode.getElementsByTagName('speciesReference');
        for i = 1:speciesReferences.getLength()
            speciesReference = speciesReferences.item(i-1);
            speciesReference.setAttribute('stoichiometry','1');
        end
    end
else
    %Create initial Node Object
    docNode = com.mathworks.xml.XMLUtils.createDocument('sbml');
    %Create Root node and define the namespace for SBML-Spatial
    docRootNode = docNode.getDocumentElement;
    docRootNode.setAttribute('xmlns','http://www.sbml.org/sbml/level3/version1');
    docRootNode.setAttribute('level', '3');%('level', '3');
    
    % D. Sullivan 9/18/13
    %probably should have an actual metaID, optional so ignoring for now
    % docRootNode.setAttribute('metaid','_000000');
    docRootNode.setAttribute('xmlns:spatial','http://www.sbml.org/sbml/level3/version1/spatial/version1');
    docRootNode.setAttribute([s,'required'],'true');
    
    %Create model node to wrap everything in
    wrapperNode = docNode.createElement('model');
    listedCompartments = 0;
end

%Create initial Node Object
%     docNode = com.mathworks.xml.XMLUtils.createDocument('sbml');
%Create Root node and define the namespace for SBML-Spatial
%     docRootNode = docNode.getDocumentElement;
docRootNode.setAttribute('xmlns','http://www.sbml.org/sbml/level3/version1/core');

docRootNode.setAttribute('level', '3');
docRootNode.setAttribute('version', '1');

%Added the 'req' namespace
%docRootNode.setAttribute('xmlns:req','http://www.sbml.org/sbml/level3/version1/requiredElements/version1');
%docRootNode.setAttribute(['req:','required'],'true');

% D. Sullivan 9/18/13
%probably should have an actual metaID, optional so ignoring for now
% docRootNode.setAttribute('metaid','_000000');
docRootNode.setAttribute('xmlns:spatial','http://www.sbml.org/sbml/level3/version1/spatial/version1');
docRootNode.setAttribute([s,'required'],'true');

wrapperNode.setAttribute('id','CellOrganizer2_7');
% wrapperNode.setAttribute('name',['CellOrganizer2_7',CSGdata.name]);
wrapperNode.setAttribute('name','CellOrganizer2_7');
wrapperNode.setAttribute('lengthUnits','um');
wrapperNode.setAttribute('areaUnits','um2');
wrapperNode.setAttribute('volumeUnits','um3');
wrapperNode.setAttribute('timeUnits','s');
wrapperNode.setAttribute('substanceUnits','molecules');


%Define units
%list
ListOfUnitDefinitions = docNode.createElement('listOfUnitDefinitions');
%substance
UnitDefinition = docNode.createElement('unitDefinition');
UnitDefinition.setAttribute('id','molecules');
ListOfUnits = docNode.createElement('listOfUnits');
unit = docNode.createElement('unit');
unit.setAttribute('kind','item');
unit.setAttribute('exponent','1');
unit.setAttribute('scale','0');
unit.setAttribute('multiplier','1');
ListOfUnits.appendChild(unit);
UnitDefinition.appendChild(ListOfUnits);
ListOfUnitDefinitions.appendChild(UnitDefinition);
%volume
UnitDefinition = docNode.createElement('unitDefinition');
UnitDefinition.setAttribute('id','um3');
ListOfUnits = docNode.createElement('listOfUnits');
unit = docNode.createElement('unit');
unit.setAttribute('kind','metre');
unit.setAttribute('exponent','3');
unit.setAttribute('scale','-6');
unit.setAttribute('multiplier','1');%'0.1');
ListOfUnits.appendChild(unit);
UnitDefinition.appendChild(ListOfUnits);
ListOfUnitDefinitions.appendChild(UnitDefinition);
%area
UnitDefinition = docNode.createElement('unitDefinition');
UnitDefinition.setAttribute('id','um2');
ListOfUnits = docNode.createElement('listOfUnits');
unit = docNode.createElement('unit');
unit.setAttribute('kind','metre');
unit.setAttribute('exponent','2');
unit.setAttribute('scale','-6');
unit.setAttribute('multiplier','1');
ListOfUnits.appendChild(unit);
UnitDefinition.appendChild(ListOfUnits);
ListOfUnitDefinitions.appendChild(UnitDefinition);
%length
UnitDefinition = docNode.createElement('unitDefinition');
UnitDefinition.setAttribute('id','um');
ListOfUnits = docNode.createElement('listOfUnits');
unit = docNode.createElement('unit');
unit.setAttribute('kind','metre');
unit.setAttribute('exponent','1');
unit.setAttribute('scale','-6');
unit.setAttribute('multiplier','1');
ListOfUnits.appendChild(unit);
UnitDefinition.appendChild(ListOfUnits);
ListOfUnitDefinitions.appendChild(UnitDefinition);
%time
UnitDefinition = docNode.createElement('unitDefinition');
UnitDefinition.setAttribute('id','s');
ListOfUnits = docNode.createElement('listOfUnits');
unit = docNode.createElement('unit');
unit.setAttribute('kind','second');
unit.setAttribute('exponent','1');
unit.setAttribute('scale','0');
unit.setAttribute('multiplier','1');
ListOfUnits.appendChild(unit);
UnitDefinition.appendChild(ListOfUnits);
ListOfUnitDefinitions.appendChild(UnitDefinition);
wrapperNode.appendChild(ListOfUnitDefinitions);


%%%Set up the Geometry node
GeowrapperNode = docNode.createElement([s,'geometry']);
GeowrapperNode.setAttribute([s,'coordinateSystem'],'cartesian');

%Define compartments
%list
inCellArray = @(given_cell_array, given_key)any(cellfun(@(x)valuesEqual(x, given_key), given_cell_array));

surfaces_compartment_pairs = cell(0, 2);
% Outside, inside
if any(strcmp(options.synthesis, {'all', 'framework', 'cell'}))
    surfaces_compartment_pairs(end+1, :) = {'EC', 'cell'};
end
if any(strcmp(options.synthesis, {'all', 'framework'}))
    surfaces_compartment_pairs(end+1, :) = {'cell', 'nucleus'};
end
for i = 1:length(models)
    model = models{i};
    model_name = model_names{i};
    model_cytonuclearflag = model_cytonuclearflags{i};
    if isfield(model, 'proteinModel')
        switch model_cytonuclearflag
            case 'cyto'
                surface_compartment_pair = {'cell', model_name};
            case 'nuc'
                surface_compartment_pair = {'nucleus', model_name};
            otherwise
                error('cytonuclearflag ''%s'' not supported', model_cytonuclearflag);
        end
        surfaces_compartment_pairs(end+1, :) = surface_compartment_pair;
    end
end
surface_list = cell(0, 1);
volumes_surfaces = containers.Map('KeyType', 'char', 'ValueType', 'char');
for i = 1:size(surfaces_compartment_pairs, 1)
    surface_compartment_pair = surfaces_compartment_pairs(i, :);
    surface_compartment_pair{1} = translateWithDefaultIdentity(name_map, surface_compartment_pair{1});
    surface_compartment_pair{2} = translateWithDefaultIdentity(name_map, surface_compartment_pair{2});
    inner_volume = surface_compartment_pair{2};
    surface_compartment_pair = sort(surface_compartment_pair);
    surfaces_compartment_pairs(i, :) = surface_compartment_pair;
    surface_compartment_pair = strjoin(surface_compartment_pair, '_');
    surface_compartment_pair = translateWithDefaultIdentity(name_map, surface_compartment_pair);
    surface_list{end+1} = surface_compartment_pair;
    volumes_surfaces(inner_volume) = surface_compartment_pair;
end
options.output.SBMLSurfaceList = surface_list;
options.output.SBMLSurfacesCompartmentPairs = surfaces_compartment_pairs;
options.output.SBMLVolumesSurfaces = volumes_surfaces;

if listedCompartments==0
    ListOfCompartments = docNode.createElement('listOfCompartments');
    compartmentlistCSG = [];
    compartmentlistMesh = [];

    if ~isempty(CSGdata)
        compartmentlistCSG = fieldnames(CSGdata);
    end
    if ~isempty(meshData)
        warning('CellOrganizer:setupSBML3','Assumes meshData(1) is frameworkMesh');
        compartmentlistMesh = arrayfun(@(x)unique(extractfield(x.list,'name')), meshData, 'UniformOutput', false);
        compartmentlistMesh = cat(2, compartmentlistMesh{:});
    end
    compartmentlist = unique([squeeze(compartmentlistCSG(:));squeeze(compartmentlistMesh(:))]);
    
    %extract compartment information
    for j = 1:length(compartmentlist)
        %actual compartment
        compartment = docNode.createElement('compartment');
        compartment.setAttribute('metaid',compartmentlist{j});
        compartment.setAttribute('id',compartmentlist{j});%CSGdata.class);
        compartment.setAttribute('name',compartmentlist{j});%CSGdata.class);
        compartment.setAttribute('spatialDimensions','3');
        compartment.setAttribute('size','50000');
        compartment.setAttribute('units','um3');%'m');
        compartment.setAttribute('constant','true');
        
        %compartment mapping
        mapping = docNode.createElement([s,'compartmentMapping']);
        mapping.setAttribute([s,'id'],[compartmentlist{j},'_compartmentMapping']);%[DomainID,CSGdata.class]);
        % Does not pass SBML validator
        if options.output.SBMLSpatialVCellCompatible
            mapping.setAttribute('id',mapping.getAttribute([s,'id']))
            mapping.setAttribute([s,'compartment'],compartmentlist{j});%[CSGdata.class]);
        end;
        mapping.setAttribute([s,'domainType'],[compartmentlist{j},'_domainType']);%DomainID);
        mapping.setAttribute([s,'unitSize'],'1');
        %add children node
        compartment.appendChild(mapping);
        
        ListOfCompartments.appendChild(compartment);
    end
    for j = 1:length(surface_list)
        %actual compartment
        compartment = docNode.createElement('compartment');
        compartment.setAttribute('id',surface_list{j});%CSGdata.class);
        compartment.setAttribute('name',surface_list{j});%CSGdata.class);
        compartment.setAttribute('spatialDimensions','2');
        compartment.setAttribute('size','50000');
        compartment.setAttribute('units','um2');%'m');
        compartment.setAttribute('constant','true');
        %compartment mapping
        mapping = docNode.createElement([s,'compartmentMapping']);
        mapping.setAttribute([s,'id'],[surface_list{j},'_compartmentMapping']);%[DomainID,CSGdata.class]);
        % Does not pass SBML validator
        if options.output.SBMLSpatialVCellCompatible
            mapping.setAttribute('id',mapping.getAttribute([s,'id']))
            mapping.setAttribute([s,'compartment'],surface_list{j});%[CSGdata.class]);
        end;
        mapping.setAttribute([s,'domainType'],[surface_list{j},'_domainType']);%DomainID);
        mapping.setAttribute([s,'unitSize'],'1');
        %add children node
        compartment.appendChild(mapping);
        
        ListOfCompartments.appendChild(compartment);
    end
    wrapperNode.appendChild(ListOfCompartments);
    
end
    

%add Species field to the files - this doesn't currently add real
%values, it is only there for demonstration. CellOrganizer doesn't
%currently support the specification of species
addspecs = 0;
if addspecs==1
    listOfSpecies = addSpecies3(docNode);
    wrapperNode.appendChild(listOfSpecies);
end

%Add listOfInitialAssignments
initAssign = 0;
if initAssign == 1
    listOfInitialAssignments = addInitAssign3(docNode);
    wrapperNode.appendChild(listOfInitialAssignments);
end
    
%%%Set up the ListOfCoordinateComponent node
ListOfCoordCompNode = docNode.createElement([s,'listOfCoordinateComponents']);
%for each dimension
dimensionNames = ['x','y','z'];
for i = 1:ndim
    CoordCompNode = docNode.createElement([s,'coordinateComponent']);
    CoordCompNode.setAttribute([s,'id'],dimensionNames(i));
    % Does not pass SBML validator
    if options.output.SBMLSpatialVCellCompatible
        CoordCompNode.setAttribute('id',CoordCompNode.getAttribute([s,'id']))
    end;
    CoordCompNode.setAttribute([s,'type'],['cartesian',upper(dimensionNames(i))]);
    CoordCompNode.setAttribute([s,'unit'],'um');
    
    %define dimensions
    minNode = docNode.createElement([s,'boundaryMin']);
    minNode.setAttribute([s,'id'],[upper(dimensionNames(i)),'min']);
    % Does not pass SBML validator
    if options.output.SBMLSpatialVCellCompatible
        minNode.setAttribute('id',minNode.getAttribute([s,'id']))
    end;
    minNode.setAttribute([s,'value'],'0');
    CoordCompNode.appendChild(minNode);
    maxNode = docNode.createElement([s,'boundaryMax']);
    maxNode.setAttribute([s,'id'],[upper(dimensionNames(i)),'max']);
    % Does not pass SBML validator
    if options.output.SBMLSpatialVCellCompatible
        maxNode.setAttribute('id',maxNode.getAttribute([s,'id']))
    end;
    switch dimensionNames(i)
        case 'x'
            maxNode.setAttribute([s,'value'], num2str(img_size_real(2)));
        case 'y'
            maxNode.setAttribute([s,'value'], num2str(img_size_real(1)));
        case 'z'
            maxNode.setAttribute([s,'value'], num2str(img_size_real(3)));
    end
    CoordCompNode.appendChild(maxNode);
    
    %add component to the list
    ListOfCoordCompNode.appendChild(CoordCompNode);
end
GeowrapperNode.appendChild(ListOfCoordCompNode);
%%%
