function [docNode,GeowrapperNode,geometryDefNodeParent,wrapperNode] = ...
    addImageObjects3(CSGdata,meshData,docNode,models,imgs,geometryDefNodeParent,GeowrapperNode,wrapperNode, options)
% ADDIMAGEOBJECTS3 Helper function that creates a SampledFieldGeometry
% instance

% Copyright (C) 2016-2019 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
% Carnegie Mellon University
% Copyright (C) 2018 Taraz Buck
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

% September 22, 2018 @tebuck Copied from addMeshObjects3, editing to 
% produce sampledFieldGeometry.

warning('off', 'CellOrganizer:addImageObjects3');
s = 'spatial:';

% Get listOfDomainTypes and listOfDomains:
ListOfDomainTypesNode = GeowrapperNode.getElementsByTagName([s,'listOfDomainTypes']).item(0);
ListOfDomainsNode = GeowrapperNode.getElementsByTagName([s,'listOfDomains']).item(0);
ListOfSampledFieldsNode = docNode.createElement([s,'listOfSampledFields']);
ListOfCoordinateComponentsNode = GeowrapperNode.getElementsByTagName([s,'listOfCoordinateComponents']).item(0);

model_names = options.output.SBMLModelNames;
model_cytonuclearflags = options.output.SBMLModelCytonuclearflags;
name_map = options.output.SBMLNameMap;
surface_list = options.output.SBMLSurfaceList;
surfaces_compartment_pairs = options.output.SBMLSurfacesCompartmentPairs;
volumes_surfaces = options.output.SBMLVolumesSurfaces;

img_size = size(imgs{1});
resolution_ijk = meshData(1).list(1).resolution([2,1,3]);
img_size_real = img_size .* resolution_ijk;

% Create a list of images for all compartments
img_list = struct('name', {}, 'ordinal', {}, 'img', {}, 'mesh', {}, 'resolution', {});
non_framework_ordinal = 3;
for i = 1:length(meshData)
    mesh_list = meshData(i).list;
    for j = 1:length(mesh_list)
        if ~isempty(mesh_list(j).img)
            img_list(end+1) = orderfields(rmfield(mesh_list(j), 'type'), img_list);
            switch(img_list(end).name)
                case translateWithDefaultIdentity(name_map, 'cell')
                    img_ordinal = 1;
                case translateWithDefaultIdentity(name_map, 'nucleus')
                    img_ordinal = 2;
                otherwise
                    img_ordinal = non_framework_ordinal;
                    non_framework_ordinal = non_framework_ordinal + 1;
            end
            img_list(end).ordinal = img_ordinal;
        end
    end
end
compartmentlistCSG = fieldnames(CSGdata);
compartmentlistCSGwithoutEC = {};
for i = 1:length(compartmentlistCSG)
    if ~strcmp(compartmentlistCSG{i}, 'EC')
        compartmentlistCSGwithoutEC{end+1} = compartmentlistCSG{i};
    end
end
warning('CellOrganizer:addImageObjects3','Do the indices correspond?');
warning('CellOrganizer:addImageObjects3','Assumes meshData(1) is frameworkMesh');
% warning('Are there better names to use?');
for i = 1:length(compartmentlistCSGwithoutEC)
    name = compartmentlistCSGwithoutEC{i};
    i2 = i + length(meshData(1).list);
    img = imgs{i2};
    if ~any(img(:))
        continue;
    end
    ordinal = length(img_list)+1;
    img_list(i2) = struct('name', name, 'ordinal', ordinal, 'img', img, 'mesh', [], 'resolution', []);
end
% Add extracellular space
ec_img = true(img_size);
for i = 1:length(imgs)
    ec_img(imgs{i} > 0) = false;
end
if any(ec_img(:))
    for i = 1:length(img_list)
        img_list(i).ordinal = img_list(i).ordinal+1;
    end
    img_list(2:end+1) = img_list;
    img_list(1) = struct('name', 'EC', 'ordinal', 1, 'img', ec_img, 'mesh', [], 'resolution', []);
end
% Sort by ordinals
[~, ordinal_order] = sort([img_list.ordinal]);
img_list = img_list(ordinal_order);

% Set size in coordinateComponents
coordinateComponents = ListOfCoordinateComponentsNode.getElementsByTagName([s,'coordinateComponent']);
for i = 1:length(coordinateComponents)
    coordinateComponent = coordinateComponents.item(i-1);
    boundaryMin = coordinateComponent.getElementsByTagName([s,'boundaryMin']).item(0);
    boundaryMax = coordinateComponent.getElementsByTagName([s,'boundaryMax']).item(0);
    boundaryMin.setAttribute([s,'value'], num2str(0));
    switch char(coordinateComponent.getAttribute([s,'type']))
        case 'cartesianX'
            boundaryMax.setAttribute([s,'value'], num2str(img_size_real(2)));
        case 'cartesianY'
            boundaryMax.setAttribute([s,'value'], num2str(img_size_real(1)));
        case 'cartesianZ'
            boundaryMax.setAttribute([s,'value'], num2str(img_size_real(3)));
    end
end


zslices = zeros(size(img_list(1).img,3),1);

use_single_sampled_field = true;

% compression = 'uncompressed';
% base64 is in draft specification 0.90 but considered invalid by both the online validator and Virtual Cell.
% compression = 'base64';
% How to encode deflated is not in the specification. Basing on https://github.com/sbmlteam/jsbml/blob/master/extensions/spatial/src/org/sbml/jsbml/ext/spatial/CompressionKind.java and https://github.com/sbmlteam/jsbml/blob/master/extensions/spatial/test/org/sbml/jsbml/xml/test/data/spatial/sampledfield_3d.xml.
% compression = 'deflated';

if options.output.SBMLSpatialUseCompression
    compression = 'deflated';
else
    compression = 'uncompressed';
end


% Create a single image for consistency

all_objects_image = [];
for i = 1:length(img_list)
    object = img_list(i);
    name = object.name;
    img = object.img;
    ordinal = object.ordinal;
    
    if isempty(all_objects_image)
        all_objects_image = uint8(img);
    end
    mask = img > 0;
    
    all_objects_image(mask) = ordinal;
end

if length(unique(all_objects_image(:)))-1 < length(img_list)
    warning('CellOrganizer:addImageObjects3','Compartments missing in indexed image');
end

all_objects_image = downsampleIntegerImage(all_objects_image, options.output.SBMLDownsampling);

%%Define new Sampled Volume object
% Sampled Field
SampledFieldNode = docNode.createElement([s,'sampledField']);
SampledFieldNodeID = ['universal','_sampledField'];
SampledFieldNode.setAttribute([s,'id'],SampledFieldNodeID);
% Does not pass SBML validator
if options.output.SBMLSpatialVCellCompatible
    SampledFieldNode.setAttribute('id',SampledFieldNode.getAttribute([s,'id']))
end;
SampledFieldNode.setAttribute([s,'numSamples1'],num2str(size(all_objects_image,2)));
SampledFieldNode.setAttribute([s,'numSamples2'],num2str(size(all_objects_image,1)));
SampledFieldNode.setAttribute([s,'numSamples3'],num2str(size(all_objects_image,3)));
% SampledFieldNode.setAttribute([s,'samplesLength'],num2str(numel(all_objects_image)));
% Validator different from specification:
% SampledFieldNode.setAttribute([s,'interpolationType'],'linear');
if options.output.SBMLSpatialVCellCompatible
    % Does not pass SBML validator
    SampledFieldNode.setAttribute([s,'interpolationType'],'nearestneighbor');
else
    SampledFieldNode.setAttribute([s,'interpolationType'],'nearestNeighbor');
end;
SampledFieldNode.setAttribute([s,'dataType'],'uint8');

[SamplesText, SamplesTextLength] = integerImageToText(all_objects_image, compression, options);
SampledFieldNode.setAttribute([s,'samplesLength'],num2str(SamplesTextLength));

% "java.lang.OutOfMemoryError: Java heap space"
SamplesTextNode = docNode.createTextNode(SamplesText);
clear SamplesText;
SampledFieldNode.setAttribute([s,'compression'], compression);
SampledFieldNode.appendChild(SamplesTextNode);

%List of Sampled Volume Objects
ListOfSampledVolumesNode = docNode.createElement([s,'listOfSampledVolumes']);

ListOfSampledFieldsNode.appendChild(SampledFieldNode);

% Split into separate images and remove blank images
img_list_consistent = img_list;
img_list_consistent_nonblank = true(size(img_list_consistent));
for i = 1:length(img_list)
    img_list_consistent(i).img_original = img_list_consistent(i).img;
    img_list_consistent(i).img = all_objects_image == img_list_consistent(i).ordinal;
    img_list_consistent_nonblank(i) = any(img_list_consistent(i).img(:));
end
img_list_consistent = img_list_consistent(img_list_consistent_nonblank);
img_list = img_list_consistent;


% Sampled Field Geometry
SampledFieldGeometryNode = docNode.createElement([s,'sampledFieldGeometry']);
SampledFieldGeometryNodeID = ['universal','_sampledFieldGeometry'];
% Validator different from specification:
SampledFieldGeometryNode.setAttribute([s,'id'],SampledFieldGeometryNodeID);
% Does not pass SBML validator
if options.output.SBMLSpatialVCellCompatible
    SampledFieldGeometryNode.setAttribute('id',SampledFieldGeometryNode.getAttribute([s,'id']))
end;
SampledFieldGeometryNode.setAttribute([s,'sampledField'],SampledFieldNodeID);
% Validator different from specification:
SampledFieldGeometryNode.setAttribute([s,'isActive'],'false');
geometryDefNodeParent.appendChild(SampledFieldGeometryNode);

warning('CellOrganizer:addImageObjects3','SBML validator (''%s'') differs from specification (''%s''), Virtual Cell agrees with specification', 'attribute "spatial:sampledValue" not allowed here', 'The required sampledValue attribute');
for i = 1:length(img_list)
    object = img_list(i);
    name = object.name;
    ordinal = object.ordinal;
    % type = object.type;
    
    % Add to listOfDomainTypes:
    DomainTypeNode = docNode.createElement([s,'domainType']);
    DomainTypeNodeID = [name,'_domainType'];
    % DomainTypeNodeID = name;
    DomainTypeNode.setAttribute([s,'id'],DomainTypeNodeID);
    % Does not pass SBML validator
    if options.output.SBMLSpatialVCellCompatible
        DomainTypeNode.setAttribute('id',DomainTypeNode.getAttribute([s,'id']))
    end;
    DomainTypeNode.setAttribute([s,'spatialDimensions'],'3');
    ListOfDomainTypesNode.appendChild(DomainTypeNode);

    % Add to listOfDomains:
    DomainNode = docNode.createElement([s,'domain']);
    % DomainNode.setAttribute([s,'id'],[name,'_domain']);
    DomainNode.setAttribute([s,'id'],name);
    % Does not pass SBML validator
    if options.output.SBMLSpatialVCellCompatible
        DomainNode.setAttribute('id',DomainNode.getAttribute([s,'id']))
    end;
    DomainNode.setAttribute([s,'domainType'],DomainTypeNodeID);
    % ListOfInteriorPoints = docNode.createElement([s,'listOfInteriorPoints']);
    % InteriorPoint = docNode.createElement([s,'interiorPoint']);
    ListOfDomainsNode.appendChild(DomainNode);
    
    ListOfInteriorPointsNode = docNode.createElement([s,'listOfInteriorPoints']);
    DomainNode.appendChild(ListOfInteriorPointsNode);
    InteriorPointNode = docNode.createElement([s,'interiorPoint']);
    ListOfInteriorPointsNode.appendChild(InteriorPointNode);
    [ipi, ipj, ipk] = find(object.img, 1, 'first');
    InteriorPointNode.setAttribute([s,'coord1'],num2str(ipj * resolution_ijk(2)));
    InteriorPointNode.setAttribute([s,'coord2'],num2str(ipi * resolution_ijk(1)));
    InteriorPointNode.setAttribute([s,'coord3'],num2str(ipk * resolution_ijk(3)));
    
    
    %assume the cell comes first
    zslicesold = zslices;
    zslices = squeeze(sum(sum(object.img,1),2));
    
    SampledVolumeNode = docNode.createElement([s,'sampledVolume']);
    SampledVolumeNode.setAttribute([s,'id'],[name,'_sampledVolume']);
    % Does not pass SBML validator
    if options.output.SBMLSpatialVCellCompatible
        SampledVolumeNode.setAttribute('id',SampledVolumeNode.getAttribute([s,'id']))
    end;
    SampledVolumeNode.setAttribute([s,'domainType'],DomainTypeNodeID);
    if use_single_sampled_field
        % Each region is stored in the single SampledFieldNode as value i on background value of zero.
        % compartment_value = i-1;
        compartment_value = ordinal;
        if options.output.SBMLSpatialVCellCompatible
            SampledVolumeNode.setAttribute([s,'sampledValue'],num2str(compartment_value));
        end
        SampledVolumeNode.setAttribute([s,'minValue'],num2str(compartment_value-0.5));
        SampledVolumeNode.setAttribute([s,'maxValue'],num2str(compartment_value+0.5));
    else
        % Each region is stored in its SampledFieldNode as value 1 on background value of zero.
        if options.output.SBMLSpatialVCellCompatible
            SampledVolumeNode.setAttribute([s,'sampledValue'],num2str(1));
        end
        SampledVolumeNode.setAttribute([s,'minValue'],num2str(0.5));
        SampledVolumeNode.setAttribute([s,'maxValue'],num2str(1.5));
    end
    
    %Get mesh points
    %remove the top and bottom slice of the nucleus fluorescence to ensure
    %that the nucleus lies within the cell.
    % if strcmpi(object.name,'NU')
    if strcmpi(object.name,'nuc')
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
    
    
    warning('CellOrganizer:addImageObjects3','Make resolution and origin apply to SampledField');
    warning('CellOrganizer:addImageObjects3','Convert addMeshObjects3 adjustments to apply to SampledField');
    
    %{
    
    %D. Sullivan 11/30/14
    %This is being replaced in favor of the iso2mesh software
    %NOTE: This adds a dependency that we may not want - we should discuss
    %FV = getMeshPoints3(object.img);
    try 
        FV = makeIso2mesh3(object.img);
    catch
        warning(['The iso2mesh package was not found.',...
            'We recommend using this package for compact high-quality meshes.',...
            'Proceeding with standard CellOrganizer method for meshing']);
        FV = getMeshPoints3(object.img);
    end
    if size(FV.faces,2)>3
        FV.faces = FV.faces(:,1:3);
    end
    
    SampledVolumeNode.setAttribute([s,'pointIndexLength'], ...
        num2str(prod(size(FV.faces))));
    
    %D. Sullivan 10/23/14
    %The mesh should be centered to fall in the range of the image.
    %Fore some reason 3.5 is being added to the z dimension.
    %This will cause the mesh to end in the wrong place.
    FV.vertices = FV.vertices.*repmat(meshData(1).resolution,size(FV.vertices,1),1);
    
    %D. Sullivan 11/10/14 - Adjust the volume and surface area for the
    %compartment
    FV.volume = sum(object.img(:)).*prod(meshData(1).resolution);
    %         FV.SA = sum(bwperim(object.img(:))).*prod(meshData(1).resolution);
    FV.SA = meshSurfaceArea(FV.vertices,FV.faces);
    
    [docNode,wrapperNode] = addVol_SurfArea3(docNode,wrapperNode,object.name,FV.volume,FV.SA);
    
    if options.debug
        save([options.temporary_results filesep ...
            'FV' object.name '_meshdata.mat'])
    end
    
    %SpatialPoints == FV.vertices
    %icaoberg COMMENT: ordinal is set to one because you only want to keep track of
    %the indices corresponding to the largest object, in this case the
    %cell. Besides the SBML schema only allows one instance of
    %SpatialPoints
    if object.ordinal==1
        SpatialPointsNode = docNode.createElement([s,'spatialPoints']);
        SpatialPointsNode.setAttribute([s,'compression'], 'uncompressed' );
        SpatialPointsNode.setAttribute([s,'arrayDataLength'], ...
            num2str(prod(size(FV.vertices))) );
        
        %icaoberg - for debugging purposes
        %vertices = docNode.createTextNode( mat2str(zeros(2,2)) );
        
        %vertices = docNode.createTextNode( mat2SBMLArrayData(FV.vertices) );
        vertices = docNode.createTextNode( [sprintf('\n'), ...
            mat2SBMLArrayData(FV.vertices, false), sprintf('\n')] );
        SpatialPointsNode.appendChild(vertices);
        SampledFieldGeometryNode.appendChild(SpatialPointsNode);
    end
    
    %Parametric Objects == FV.faces
    %icaoberg - for debugging purposes
    %faces = docNode.createTextNode( mat2SBMLArrayData(FV.faces) );
    faces = docNode.createTextNode( [sprintf('\n'), ...
        mat2SBMLArrayData(FV.faces, false), sprintf('\n')] );
    
    %}
    
    
    if ~use_single_sampled_field
        object_image = object.img;
        object_image = downsampleIntegerImage(object_image, options.output.SBMLDownsampling);
        
        % Sampled Field
        SampledFieldNode = docNode.createElement([s,'sampledField']);
        SampledFieldNodeID = [name,'_sampledField'];
        SampledFieldNode.setAttribute([s,'id'],SampledFieldNodeID);
        % Does not pass SBML validator
        if options.output.SBMLSpatialVCellCompatible
            SampledFieldNode.setAttribute('id',SampledFieldNode.getAttribute([s,'id']))
        end;
        SampledFieldNode.setAttribute([s,'numSamples1'],num2str(size(object_image,2)));
        SampledFieldNode.setAttribute([s,'numSamples2'],num2str(size(object_image,1)));
        SampledFieldNode.setAttribute([s,'numSamples3'],num2str(size(object_image,3)));
        % Validator different from specification:
        % SampledFieldNode.setAttribute([s,'interpolationType'],'linear');
        if options.output.SBMLSpatialVCellCompatible
            % Does not pass SBML validator
            SampledFieldNode.setAttribute([s,'interpolationType'],'nearestneighbor');
        else
            SampledFieldNode.setAttribute([s,'interpolationType'],'nearestNeighbor');
        end;
        SampledFieldNode.setAttribute([s,'dataType'],'uint8');
        
        %List of Sampled Volume Objects
        ListOfSampledVolumesNode = docNode.createElement([s,'listOfSampledVolumes']);
        
        [SamplesText, SamplesTextLength] = integerImageToText(object_image, compression, options);
        SampledFieldNode.setAttribute([s,'samplesLength'],num2str(SamplesTextLength));
        
        % "java.lang.OutOfMemoryError: Java heap space"
        SamplesTextNode = docNode.createTextNode(SamplesText);
        clear SamplesText;
        SampledFieldNode.setAttribute([s,'compression'], compression);
        SampledFieldNode.appendChild(SamplesTextNode);
        
        ListOfSampledFieldsNode.appendChild(SampledFieldNode);
    end
    
    ListOfSampledVolumesNode.appendChild(SampledVolumeNode);
    SampledFieldGeometryNode.appendChild(ListOfSampledVolumesNode);
    % geometryDefNode.appendChild(SampledFieldGeometryNode);
end

%wrapperNode.appendChild(ListOfCompartments);
GeowrapperNode.appendChild(ListOfSampledFieldsNode);



function Out = fBase64_enc(In)
% From DataHash.m:
% Author: Jan Simon, Heidelberg, (C) 2011 matlab.THISYEAR(a)nMINUSsimon.de
% License: BSD (use/copy/change/redistribute on own risk, mention the author)
% Encode numeric vector of UINT8 values to base64 string.

Pool = [65:90, 97:122, 48:57, 43, 47];  % [0:9, a:z, A:Z, +, /]
v8   = [128; 64; 32; 16; 8; 4; 2; 1];
v6   = [32, 16, 8, 4, 2, 1];

In  = reshape(In, 1, []);
X   = rem(floor(In(ones(8, 1), :) ./ v8(:, ones(length(In), 1))), 2);
Y   = reshape([X(:); zeros(6 - rem(numel(X), 6), 1)], 6, []);
Out = char(Pool(1 + v6 * Y));



function wrappedString = wrapString(stringData, wrapWidth, newlineString, separatorRegexp)
% tebuck

if nargin < 2
    wrapWidth = 72;
end
if nargin < 3
    newlineString = sprintf('\n');
end
if nargin < 4
    separatorRegexp = '\s+';
end

if strcmp(separatorRegexp, '')
    wrappedString = stringData;
    wrappedStringNumRows = ceil(length(wrappedString) / wrapWidth);
    wrappedString = [mat2cell(reshape(wrappedString(1:(wrappedStringNumRows - 1) * wrapWidth), wrapWidth, [])', ones(wrappedStringNumRows - 1, 1)); {wrappedString((wrappedStringNumRows - 1) * wrapWidth:end)}];
    wrappedString = strcat(wrappedString, {newlineString});
    wrappedString = [wrappedString{:}];
else
    wrappedStringCell = {''};
    [separatorStartPositions, separatorEndPositions] = regexp(stringData, separatorRegexp, 'start', 'end');
    separatorLengths = separatorEndPositions-separatorStartPositions+1;
    currentPosition = 1;
    currentSeparatorPositionIndex = 0;
    while currentPosition <= length(stringData)
        currentEndPosition = currentPosition;
        nextPosition = currentPosition;
        while true
            if currentSeparatorPositionIndex == length(separatorStartPositions)
                currentEndPosition = length(stringData)+1;
                nextPosition = currentEndPosition+1;
                break;
            end
            nextWrapStartPosition = separatorStartPositions(currentSeparatorPositionIndex+1);
            nextWrapEndPosition = separatorEndPositions(currentSeparatorPositionIndex+1);
            nextWrapSkip = separatorLengths(currentSeparatorPositionIndex+1);
            if nextWrapStartPosition - currentPosition < wrapWidth || currentEndPosition == currentPosition
                currentEndPosition = nextWrapStartPosition;
                nextPosition = currentEndPosition+nextWrapSkip;
                currentSeparatorPositionIndex = currentSeparatorPositionIndex+1;
            else
                break;
            end
        end
        wrappedStringCell{end+1} = stringData(currentPosition:currentEndPosition-1);
        currentPosition = nextPosition;
    end
    wrappedString = strjoin(wrappedStringCell, newlineString);
end



function img = downsampleIntegerImage(img, SBMLDownsampling)

if any(SBMLDownsampling < 1)
    classfunc = str2func(class(img));
    % img = classfunc(round(imresize(single(img), SBMLDownsampling)));
    img = classfunc(tp_imresize3d(img,round(size(img).*SBMLDownsampling),'nearest'));
end



function [CompressedSamplesText, CompressedSamplesTextLength] = integerImageToText(img, compression, options)

SamplesBytes = img;
SamplesBytes = permute(SamplesBytes, [2, 1, 3]);
SamplesBytes = SamplesBytes(:);

if strcmp(compression, 'uncompressed')
    SamplesText = mat2SBMLArrayData(uint8(SamplesBytes), true);
    SamplesTextLength = numel(SamplesBytes);

    CompressedSamplesText = SamplesText;
    CompressedSamplesTextLength = SamplesTextLength;
    
elseif strcmp(compression, 'base64')
    error('Invalid value for compression ''%s'' (in specifications but not supported by validator or Virtual Cell)', compression);
    
    SamplesText = mat2SBMLArrayData(uint8(SamplesBytes), true);
    SamplesTextLength = numel(SamplesBytes);

    % base64encode not available until R2016b
    % SamplesText = base64encode(all_objects_image(:));
    % SamplesText = fBase64_enc(uint8(all_objects_image(:)));
    CompressedSamplesText = fBase64_enc(SamplesText);
    CompressedSamplesText = wrapString(CompressedSamplesText);
    error('Compression ''%s'' not implemented (should inputStream be ASCII as with uncompressed?)', compression);
    CompressedSamplesText = [sprintf('\n'), CompressedSamplesText, sprintf('\n')];
    CompressedSamplesTextLength = SamplesTextLength;
    
elseif strcmp(compression, 'deflated')
    [CompressedSamplesText, CompressedSamplesTextLength] = SBMLDeflate(SamplesBytes, options);
    
    %{
    % Test
    [DecompressedSamplesBytes] = SBMLInflate(CompressedSamplesText, options);
    CompressedSamplesTextFirst100Characters = CompressedSamplesText(1:100)
    SamplesBytesFirst100 = SamplesBytes(1:100)'
    DecompressedSamplesBytesFirst100 = DecompressedSamplesBytes(1:100)
    matches = all(SamplesBytesFirst100 == DecompressedSamplesBytesFirst100)
    mean_match = mean(SamplesBytesFirst100 == DecompressedSamplesBytesFirst100)
    error('DEBUG');
    %}
    
else
    error('Invalid value for compression ''%s''', compression);
end



function [CompressedSamplesText, CompressedSamplesTextLength] = SBMLDeflate(SamplesBytes, options)

format = 'deflate';

SamplesBytes = uint8(SamplesBytes);

% https://undocumentedmatlab.com/blog/savezip-utility
if options.output.SBMLSpatialVCellCompatible
    inputStream = uint8(SamplesBytes);
else
    SamplesText = strjoin(arrayfun(@int2str, uint8(SamplesBytes), 'UniformOutput', false));
    inputStream = uint8(SamplesText);
end
outputStream = java.io.ByteArrayOutputStream();
switch format
    case 'deflate'
        % Does not work in VCell
        compressedInputStream = java.util.zip.DeflaterOutputStream(outputStream);
        % Does not work in VCell
        % compressedInputStream = java.util.zip.DeflaterOutputStream(outputStream, java.util.zip.Deflater(java.util.zip.Deflater.DEFAULT_COMPRESSION, true));
    case 'gzip'
        compressedInputStream = java.util.zip.GZIPOutputStream(outputStream);
    otherwise
        error('Unknown format ''%s''', format);
end
compressedInputStream.write(inputStream, 0, numel(inputStream));
compressedInputStream.finish();
compressedBytes = outputStream.toByteArray();
% This does not work with Virtual Cell
compressedBytes = typecast(compressedBytes, 'uint8');
compressedInputStream.close();
outputStream.close();
clear compressedInputStream;
clear outputStream;
clear inputStream;
CompressedSamplesText = strjoin(arrayfun(@int2str, compressedBytes, 'UniformOutput', false));
CompressedSamplesTextLength = length(compressedBytes);
clear compressedBytes;



function [DecompressedSamplesBytes] = SBMLInflate(CompressedSamplesText, options)

format = 'deflate';

compressedBytes = str2double(strsplit(CompressedSamplesText));
% This does not work with Virtual Cell
% compressedBytes = uint8(compressedBytes);
inputStream = java.io.ByteArrayInputStream(compressedBytes);
switch format
    case 'deflate'
        % Does not work in VCell
        decompressedInputStream = java.util.zip.InflaterInputStream(inputStream);
        % Does not work in VCell
        % decompressedInputStream = java.util.zip.InflaterInputStream(inputStream, java.util.zip.Deflater(java.util.zip.Deflater.DEFAULT_COMPRESSION, true));
    case 'gzip'
        decompressedInputStream = java.util.zip.GZIPInputStream(inputStream);
    otherwise
        error('Unknown format ''%s''', format);
end
DecompressedSamplesBytesText = org.apache.commons.io.IOUtils.toByteArray(decompressedInputStream)';
decompressedInputStream.close();
inputStream.close();
clear decompressedInputStream;
clear inputStream;
clear compressedBytes;
DecompressedSamplesBytesTextFirst100 = DecompressedSamplesBytesText(1:100)
DecompressedSamplesBytesText = typecast(DecompressedSamplesBytesText, 'uint8');
if options.output.SBMLSpatialVCellCompatible
    DecompressedSamplesBytes = DecompressedSamplesBytesText;
else
    DecompressedSamplesBytesText = char(DecompressedSamplesBytesText);
    DecompressedSamplesBytes = uint8(str2double(strsplit(DecompressedSamplesBytesText)));
end

