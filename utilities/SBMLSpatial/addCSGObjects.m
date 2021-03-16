function [docNode,GeowrapperNode,geometryDefNode,wrapperNode] = addCSGObjects(CSGdata,docNode,geometryDefNode,GeowrapperNode,wrapperNode,options)

s = 'spatial:';
% listedCompartments = 0;

rotation_prefix = 'rotate';
if options.output.SBMLSpatialVCellCompatible
    % rotateAxis* is currently invalid SBML, but VCell still needs it
    rotation_prefix = 'rotateAxis';
end

%%%Set up the ListOfDomainTypes node
ListOfDomainTypesNode = docNode.createElement([s,'listOfDomainTypes']);
%get list of domain types from CSG data list
%     domainlist = unique(extractfield(CSGdata.list,'name'));
domainlist = fieldnames(CSGdata);
domainlist2 = {};
for i = 1:length(domainlist)
    if ~strcmp(domainlist{i}, 'primitiveOnly')
        domainlist2{end+1} = domainlist{i};
    end
end
domainlist = domainlist2;
for j = 1:length(domainlist)
    DomainTypeNode = docNode.createElement([s,'domainType']);
    DomainTypeNode.setAttribute([s,'id'],domainlist{j});
    % Does not pass SBML validator
    if options.output.SBMLSpatialVCellCompatible
        DomainTypeNode.setAttribute('id',DomainTypeNode.getAttribute([s,'id']))
    end
    DomainTypeNode.setAttribute([s,'spatialDimensions'],'3');
    ListOfDomainTypesNode.appendChild(DomainTypeNode);
    
    GeowrapperNode.appendChild(ListOfDomainTypesNode);
    %%%Set up ListOfDomains
    ListOfDomains = docNode.createElement([s,'listOfDomains']);
    for i=1:length(CSGdata.(domainlist{j}).list)
        %     if (j~=67&&j~=length(CSGdata.list))
        %         continue
        %     end
        object = CSGdata.(domainlist{j}).list(i);
        name = object.name;
        type = object.type;
        DomainNode = docNode.createElement([s,'domain']);
        DomainNode.setAttribute([s,'id'],[name,num2str(i-1)]);%[DomainID,num2str(j-1)])%'0']);
        % Does not pass SBML validator
        if options.output.SBMLSpatialVCellCompatible
            DomainNode.setAttribute('id',DomainNode.getAttribute([s,'id']))
        end
        DomainNode.setAttribute([s,'domainType'],name);
        ListOfInteriorPoints = docNode.createElement([s,'listOfInteriorPoints']);
        InteriorPoint = docNode.createElement([s,'interiorPoint']);
        if isfield(object,'position')
            position = object.position;
        else
            position = [0,0,0];
        end
        InteriorPoint.setAttribute([s,'coord1'],num2str(position(1)));%'0');
        InteriorPoint.setAttribute([s,'coord2'],num2str(position(2)));%'0');
        InteriorPoint.setAttribute([s,'coord3'],num2str(position(3)));%'0');
        ListOfInteriorPoints.appendChild(InteriorPoint);
        DomainNode.appendChild(ListOfInteriorPoints);
        ListOfDomains.appendChild(DomainNode);
    end
%     ordinals(j) = CSGdata.(domainlist{j}).list.ordinal;
end
GeowrapperNode.appendChild(ListOfDomains);

%%%Define all the CSG (Primitive) type geometries
CSGeometryNode = docNode.createElement([s,'csGeometry']);
CSGeometryNode.setAttribute([s,'id'], 'CSG_Geometry1');
% Does not pass SBML validator
if options.output.SBMLSpatialVCellCompatible
    CSGeometryNode.setAttribute('id',CSGeometryNode.getAttribute([s,'id']))
end
CSGeometryNode.setAttribute([s,'isActive'],'true');
% geometryDefNode.appendChild(CSGeometryNode);

ListOfCSGObjectsNode = docNode.createElement([s,'listOfCSGObjects']);

bounding_box_object = [];
for i = 1:length(CSGdata.EC.list)
    csg_object = CSGdata.EC.list(i);
    found_bounding_box = strcmp(csg_object.name, 'EC') && strcmp(csg_object.type, 'cube');
    if ~found_bounding_box
        continue;
    end
    bounding_box_object = csg_object;
    if found_bounding_box
        break;
    end
end

%add each CSGObject
for k = 1:length(domainlist)
    for j=1:length(CSGdata.(domainlist{k}).list)
        %     if (j~=67&&j~=length(CSGdata.list))
        %         continue
        %     end
        object = CSGdata.(domainlist{k}).list(j);
        name = object.name;
        type = object.type;
        object.ordinal = CSGdata.(domainlist{k}).ordinal;
        
        CSGObjectNode = docNode.createElement([s,'csgObject']);
        %     CSGObjectNode.setAttribute([s,'spatialID'],['Sp_', name]);
        CSGObjectNode.setAttribute([s,'id'],name);
        % Does not pass SBML validator
        if options.output.SBMLSpatialVCellCompatible
            CSGObjectNode.setAttribute('id',CSGObjectNode.getAttribute([s,'id']))
        end
        CSGObjectNode.setAttribute([s,'domainType'],name);
        
        if isfield(object, 'ordinal')
            ordinal = num2str(object.ordinal);
            CSGObjectNode.setAttribute([s,'ordinal'], ordinal);
          
          %%This is now done in the getBox.m method as it must be consistent with the volume parameters saved  
        end
        
        CSGTransformationNode = docNode.createElement([s,'csgHomogeneousTransformation']);
        
        if isfield(object, 'position')
            position = object.position;
            if options.output.SBMLFlipXToAlign
                % Flip to align with image output
                % position(1) = bounding_box_object.position(1) - (position(1) - bounding_box_object.position(1));
                max_size_voxels = size(options.cell);
                max_size_voxels_xyz = max_size_voxels([2, 1, 3]);
                max_size = max_size_voxels_xyz .* options.resolution.cubic;
                position(1) = max_size(1) - (position(1) - 1 * options.resolution.cubic(1));
            end
            CSGTranslationNode = docNode.createElement([s,'csgTranslation']);
            CSGTranslationNode.setAttribute([s,'id'],'translation');
            % Does not pass SBML validator
            if options.output.SBMLSpatialVCellCompatible
                CSGTranslationNode.setAttribute('id',CSGTranslationNode.getAttribute([s,'id']))
            end
            CSGTranslationNode.setAttribute([s,'translateX'], num2str(position(1)));
            CSGTranslationNode.setAttribute([s,'translateY'], num2str(position(2)));
            CSGTranslationNode.setAttribute([s,'translateZ'], num2str(position(3)));
        else
            error('No position specified. Unable to create SBML Spatial file.');
        end
        
        if isfield(object, 'rotation')
            if ~options.output.SBMLFlipXToAlign
                warning('Object rotations are correct for options.output.SBMLFlipXToAlign = true only');
            end
            rotation = object.rotation;
            rotation_matrix = object.rotationmatrix;
            CSGRotationNode = docNode.createElement([s,'csgRotation']);
            CSGRotationNode.setAttribute([s,'id'], 'rotation');
            % Does not pass SBML validator
            if options.output.SBMLSpatialVCellCompatible
                CSGRotationNode.setAttribute('id',CSGRotationNode.getAttribute([s,'id']))
                CSGRotationNode.setAttribute([s,rotation_prefix,'X'], num2str(deg2rad(rotation(1))));
                CSGRotationNode.setAttribute([s,rotation_prefix,'Y'], num2str(deg2rad(rotation(2))));
                CSGRotationNode.setAttribute([s,rotation_prefix,'Z'], num2str(deg2rad(rotation(3))));
            end
            % rotation is in degrees (objrotvec in tp_gengaussobjimg.m)
            % rotation_matrix = eulerAnglesToRotation3d(rotation);
            % Rotate to align with image output
            rotation_matrix = createRotationOz(pi/2) * rotation_matrix;
            [rotation_axis, rotation_angle] = rotation3dAxisAndAngle(rotation_matrix);
            rotation_axis = rotation_axis(4:6);
            CSGRotationNode.setAttribute([s,rotation_prefix,'X'], num2str(rotation_axis(1)));
            CSGRotationNode.setAttribute([s,rotation_prefix,'Y'], num2str(rotation_axis(2)));
            CSGRotationNode.setAttribute([s,rotation_prefix,'Z'], num2str(rotation_axis(3)));
            CSGRotationNode.setAttribute([s,'rotateAngleInRadians'], num2str(rotation_angle));
        else
            error('No rotation specified. Unable to create SBML Spatial file.');
        end
        
        if isfield(object, 'scale')
            objsize = object.scale;%.*resolution;
            CSGScaleNode = docNode.createElement([s,'csgScale']);
            CSGScaleNode.setAttribute([s,'id'], 'scale');
            % Does not pass SBML validator
            if options.output.SBMLSpatialVCellCompatible
                CSGScaleNode.setAttribute('id',CSGScaleNode.getAttribute([s,'id']))
            end
            CSGScaleNode.setAttribute([s,'scaleX'], num2str(objsize(1)));
            CSGScaleNode.setAttribute([s,'scaleY'], num2str(objsize(2)));
            CSGScaleNode.setAttribute([s,'scaleZ'], num2str(objsize(3)));
        else
            error('No scale specified. Unable to create SBML Spatial file.');
        end
        
        CSGPrimitiveNode = docNode.createElement([s,'csgPrimitive']);
        CSGPrimitiveNode.setAttribute([s,'id'],type);
        % Does not pass SBML validator
        if options.output.SBMLSpatialVCellCompatible
            CSGPrimitiveNode.setAttribute('id',CSGPrimitiveNode.getAttribute([s,'id']))
        end
        % CSGPrimitiveNode.setAttribute([s,'spatialId'],'cube');
        if strcmpi('sphere',type)
            CSGPrimitiveNode.setAttribute([s,'primitiveType'],'sphere');
        elseif strcmpi('cube',type)
            CSGPrimitiveNode.setAttribute([s,'primitiveType'],'cube');
        else
            warning('Unrecognized primitive type specified, defaulting to cube');
            CSGPrimitiveNode.setAttribute([s,'primitiveType'],'cube');
        end
        
        CSGTransformationNode.appendChild(CSGRotationNode);
        CSGTransformationNode.appendChild(CSGTranslationNode);
        CSGTransformationNode.appendChild(CSGScaleNode);
        CSGObjectNode.appendChild(CSGTransformationNode);
        CSGObjectNode.appendChild(CSGPrimitiveNode);
        ListOfCSGObjectsNode.appendChild(CSGObjectNode);
        
    end
end

if ListOfCSGObjectsNode.getChildNodes().getLength() > 0
    CSGeometryNode.appendChild(ListOfCSGObjectsNode);
end
if CSGeometryNode.getChildNodes().getLength() > 0
    geometryDefNode.appendChild(CSGeometryNode);
end
