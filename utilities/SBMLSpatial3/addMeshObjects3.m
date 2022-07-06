function [docNode,GeowrapperNode,geometryDefNode,wrapperNode] = ...
    addMeshObjects3(meshData,docNode,geometryDefNode,GeowrapperNode,wrapperNode, options)
% ADDMESHOBJECTS3 Helper function that creates a ParametricGeometry
% instance

% Copyright (C) 2016 Murphy Lab
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

% November 29, 2016 @icaoberg Updated method to save number of items in 
%        spatialPoints and parametricGeometry instances


debug = false;
% debug = true;

if ~debug
    warning('off', 'CellOrganizer:addMeshObjects3');
end


s = 'spatial:';

surface_list = options.output.SBMLSurfaceList;
surfaces_compartment_pairs = options.output.SBMLSurfacesCompartmentPairs;
volumes_surfaces = options.output.SBMLVolumesSurfaces;

% Get listOfDomainTypes and listOfDomains:
ListOfDomainTypesNode = GeowrapperNode.getElementsByTagName([s,'listOfDomainTypes']).item(0);
ListOfDomainsNode = GeowrapperNode.getElementsByTagName([s,'listOfDomains']).item(0);

%%Define all the Mesh objects
warning('CellOrganizer:addMeshObjects3','Assumes meshData(1) is frameworkMesh');
zslices = zeros(size(meshData(1).list(1).img,3),1);

all_vertices = zeros(0, 3);

k = 1;
for i = 1:length(meshData)
    for j = 1:length(meshData(i).list)
        % Parametric Geometry
        ParaGeometryNode = docNode.createElement([s,'parametricGeometry']);
        ParaGeometryNode.setAttribute([s,'id'],sprintf('parametricGeometry%d',k));
        % Does not pass SBML validator
        if options.output.SBMLSpatialVCellCompatible
            ParaGeometryNode.setAttribute('id',ParaGeometryNode.getAttribute([s,'id']))
        end;
        ParaGeometryNode.setAttribute([s,'isActive'],'false');

        % geometryDefNode.appendChild(ParaGeometryNode);

        object = meshData(i).list(j);
        name = object.name;
        type = object.type;
        object_mesh = object.mesh;
        object_resolution = object.resolution;
        
        % Add to listOfDomainTypes:
        ListOfDomainTypesNode = GeowrapperNode.getElementsByTagName([s,'listOfDomainTypes']).item(0);
        DomainTypeNode = docNode.createElement([s,'domainType']);
        DomainTypeNode.setAttribute([s,'id'],[name,'_domainType']);
        % Does not pass SBML validator
        if options.output.SBMLSpatialVCellCompatible
            DomainTypeNode.setAttribute('id',DomainTypeNode.getAttribute([s,'id']))
        end;
        DomainTypeNode.setAttribute([s,'spatialDimensions'],'3');
        ListOfDomainTypesNode.appendChild(DomainTypeNode);

        % Add to listOfDomains:
        DomainNode = docNode.createElement([s,'domain']);
        DomainNode.setAttribute([s,'id'],name);
        % Does not pass SBML validator
        if options.output.SBMLSpatialVCellCompatible
            DomainNode.setAttribute('id',DomainNode.getAttribute([s,'id']))
        end;
        DomainNode.setAttribute([s,'domainType'],[name,'_domainType']);
        ListOfInteriorPoints = docNode.createElement([s,'listOfInteriorPoints']);
        InteriorPoint = docNode.createElement([s,'interiorPoint']);
        position = object_mesh.vertices(1, :);
        position_found = false;
        for index1 = 1:size(object_mesh.vertices, 1)
            for index2 = 1:size(object_mesh.vertices, 1)
                position = mean(object_mesh.vertices([index1, index2], :), 1);
                position_found = meshContainsPoint(object_mesh, position);
                if position_found
                    break;
                end
            end
            if position_found
                break;
            end
        end
        if ~position_found
            error('Could not find interior point of mesh');
        end
        InteriorPoint.setAttribute([s,'coord1'],num2str(position(1)));
        InteriorPoint.setAttribute([s,'coord2'],num2str(position(2)));
        InteriorPoint.setAttribute([s,'coord3'],num2str(position(3)));
        ListOfInteriorPoints.appendChild(InteriorPoint);
        DomainNode.appendChild(ListOfInteriorPoints);
        ListOfDomainsNode.appendChild(DomainNode);
        
        %assume the cell comes first
        zslicesold = zslices;
        zslices = squeeze(sum(sum(object.img,1),2));
        
        ParaObjectNode = docNode.createElement([s,'parametricObject']);
        % ParaObjectNode.setAttribute([s,'id'],[name,num2str(1),'_mesh']);%[name]);%['Sp_', name]);
        ParaObjectNode.setAttribute([s,'id'],name);%[name]);%['Sp_', name]);
        % Does not pass SBML validator
        if options.output.SBMLSpatialVCellCompatible
            ParaObjectNode.setAttribute('id',ParaObjectNode.getAttribute([s,'id']))
        end;
        %for now we will assume we are only dealing with triangulated meshes
        ParaObjectNode.setAttribute([s,'polygonType'],'triangle');
        %D.Sullivan 4/15/14 - fixed naming bug.
        ParaObjectNode.setAttribute([s,'domainType'],[name,'_domainType']);
        ParaObjectNode.setAttribute([s,'compression'],'uncompressed');
        
        %Get mesh points
        %remove the top and bottom slice of the nucleus fluorescence to ensure
        %that the nucleus lies within the cell.
        if strcmpi(object.name,'NU')
            slices = find(max(max(object.img))>0);
            botslice = min(slices);
            topslice = max(slices);
            %D. Sullivan 12/2/14 - this actually doesn't do anything. need to
            %remove the top and bottom slices of nuc. 
    %         object.img(:,:,1) = object.img(:,:,botslice).*0;
            object.img(:,:,botslice) = object.img(:,:,botslice).*0;
    %         object.img(:,:,end) = object.img(:,:,topslice).*0;
            object.img(:,:,topslice) = object.img(:,:,topslice).*0;
            
        end
        
        image_to_contour = object.img;
        if options.output.SBMLFlipXToAlign
            % Flip in X to align with image output
            image_to_contour = flipdim(image_to_contour, 2);
        end
        
        use_object_mesh = isstruct(object_mesh) && options.output.SBMLSpatialUseAnalyticMeshes;
        if use_object_mesh
            FV = object_mesh;
            FV.vertices = FV.vertices ./ repmat(object_resolution, size(FV.vertices, 1), 1);
            %{
            if options.output.SBMLFlipXToAlign
                % Flip in X to align with image output
                % FV.vertices(:, 1) = size(object.img, 2) + 1 - FV.vertices(:, 1);
                max_size_voxels = size(options.cell);
                max_size_voxels_xyz = max_size_voxels([2, 1, 3]);
                max_size = max_size_voxels_xyz .* options.resolution.cubic;
                FV.vertices(:, 1) = max_size_voxels_xyz(1) - (FV.vertices(:, 1) - 1);
                FV.faces(:, 2:3) = FV.faces(:, [3, 2]);
            end
            %}
        else
            %D. Sullivan 11/30/14
            %This is being replaced in favor of the iso2mesh software
            %NOTE: This adds a dependency that we may not want - we should discuss
            %FV = getMeshPoints3(image_to_contour);
            try 
                % warning('makeIso2mesh3 temporarily disabled because it is not working by default'); error();
                FV = makeIso2mesh3(image_to_contour);
            catch
                warning(['The iso2mesh package was not found.',...
                    'We recommend using this package for compact high-quality meshes.',...
                    'Proceeding with standard CellOrganizer method for meshing']);
                if isempty(image_to_contour)
                    FV = struct('vertices', zeros(0, 3), 'faces', zeros(0, 3));
                else
                    FV = getMeshPoints3(image_to_contour, [], options.output.SBMLDownsampling, options);
                end
            end
            if size(FV.faces,2) == 4
                FV.faces = [FV.faces(:, 1:3); FV.faces(:, [1, 3, 4])];
            end
        end
        
        if size(FV.vertices,1) == 0
            continue
        end
        
        ParaObjectNode.setAttribute([s,'pointIndexLength'], ...
            num2str(prod(size(FV.faces))));
        
        %D. Sullivan 10/23/14
        %The mesh should be centered to fall in the range of the image.
        %Fore some reason 3.5 is being added to the z dimension.
        %This will cause the mesh to end in the wrong place.
        FV.vertices = FV.vertices.*repmat(object_resolution,size(FV.vertices,1),1);
        
        %D. Sullivan 11/10/14 - Adjust the volume and surface area for the
        %compartment
        FV.volume = sum(object.img(:)).*prod(object_resolution);
        %         FV.SA = sum(bwperim(object.img(:))).*prod(object_resolution);
        FV.SA = meshSurfaceArea(FV.vertices,FV.faces);
        
        [docNode,wrapperNode] = addVol_SurfArea3(docNode,wrapperNode,object.name,FV.volume,FV.SA);
        
        if options.debug
            save([options.temporary_results filesep ...
                'FV' object.name '_meshdata.mat'])
        end
        
        % Use annotation to specify surface class
        if volumes_surfaces.isKey(name)
            object_annotation = docNode.createElement(['annotation']);
            cb = 'cellblender:';
            object_annotation.setAttribute(['xmlns:','cellblender'],'http://mcell.org/CellBlender');
            object_annotation_gsc = docNode.createElement([cb,'globalSurfaceClass']);
            object_annotation_gsc.setAttribute([cb,'name'],volumes_surfaces(name));
            object_annotation.appendChild(object_annotation_gsc);
            ParaObjectNode.appendChild(object_annotation);
        end
        

        %List of Parametric Objects
        ListOfParaObjectsNode = docNode.createElement([s,'listOfParametricObjects']);
        ListOfParaObjectsNode.appendChild(ParaObjectNode);
        ParaGeometryNode.appendChild(ListOfParaObjectsNode);
        
        %SpatialPoints
        SpatialPointsNode = docNode.createElement([s,'spatialPoints']);
        SpatialPointsNode.setAttribute([s,'compression'], 'uncompressed' );
        SpatialPointsNode.setAttribute([s,'arrayDataLength'], ...
            num2str(prod(size(FV.vertices))) );
        vertices = docNode.createTextNode( [sprintf('\n'), ...
            mat2SBMLArrayData(FV.vertices, false), sprintf('\n')] );
        SpatialPointsNode.appendChild(vertices);
        ParaGeometryNode.appendChild(SpatialPointsNode);

        %wrapperNode.appendChild(ListOfCompartments);
        geometryDefNode.appendChild(ParaGeometryNode);
        
        k = k + 1;
    end
end

