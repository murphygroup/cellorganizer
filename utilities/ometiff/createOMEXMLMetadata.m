function metadata = createOMEXMLMetadata(I, parameters, varargin)
% CREATEMINIMALOMEXMLMETADATA Create an OME-XML metadata object from an input matrix
%
%    createMinimalOMEXMLMetadata(I) creates an OME-XML metadata object from
%    an input 5-D array. Minimal metadata information is stored such as the
%    pixels dimensions, dimension order and type. The output object is a
%    metadata object of type loci.formats.ome.OMEXMLMetadata.
%
%    createMinimalOMEXMLMetadata(I, dimensionOrder) specifies the dimension
%    order of the input matrix. Default valuse is XYZCT.
%
%    Examples:
%
%        metadata = createMinimalOMEXMLMetadata(zeros(100, 100));
%        metadata = createMinimalOMEXMLMetadata(zeros(10, 10, 2), 'XYTZC');
%
% See also: BFSAVE

% OME Bio-Formats package for reading and converting biological file formats.
%
% Copyright (C) 2012 - 2016 Open Microscopy Environment:
%   - Board of Regents of the University of Wisconsin-Madison
%   - Glencoe Software, Inc.
%   - University of Dundee
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as
% published by the Free Software Foundation, either version 2 of the
% License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc.,
% 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

% Not using the inputParser for first argument as it copies data

% additional author: Xin Lu (xlu2@andrew.cmu.edu)
assert(isnumeric(I), 'First argument must be numeric');

% Input check
ip = inputParser;
ip.addOptional('dimensionOrder', 'XYZCT', @(x) ismember(x, getDimensionOrders()));
ip.parse(varargin{:});

% Create metadata
toInt = @(x) javaObject('ome.xml.model.primitives.PositiveInteger', ...
    javaObject('java.lang.Integer', x));
OMEXMLService = javaObject('loci.formats.services.OMEXMLServiceImpl');
metadata = OMEXMLService.createOMEXMLMetadata();
metadata.createRoot();
metadata.setImageID('Image:0', 0);
metadata.setPixelsID('Pixels:0', 0);

if is_octave()
    java_true = java_get('java.lang.Boolean', 'TRUE');
else
    java_true = java.lang.Boolean.TRUE;
end
metadata.setPixelsBigEndian(java_true, 0);

% Set dimension order
dimensionOrderEnumHandler = javaObject('ome.xml.model.enums.handlers.DimensionOrderEnumHandler');
dimensionOrder = dimensionOrderEnumHandler.getEnumeration(ip.Results.dimensionOrder);
metadata.setPixelsDimensionOrder(dimensionOrder, 0);

% Set pixels type
pixelTypeEnumHandler = javaObject('ome.xml.model.enums.handlers.PixelTypeEnumHandler');
if strcmp(class(I), 'single')
    pixelsType = pixelTypeEnumHandler.getEnumeration('float');
else
    pixelsType = pixelTypeEnumHandler.getEnumeration(class(I));
end
metadata.setPixelsType(pixelsType, 0);

% Read pixels size from image and set it to the metadat
sizeX = size(I, 2);
sizeY = size(I, 1);
sizeZ = size(I, find(ip.Results.dimensionOrder == 'Z'));
sizeC = size(I, find(ip.Results.dimensionOrder == 'C'));
sizeT = size(I, find(ip.Results.dimensionOrder == 'T'));
metadata.setPixelsSizeX(toInt(sizeX), 0);
metadata.setPixelsSizeY(toInt(sizeY), 0);
metadata.setPixelsSizeZ(toInt(sizeZ), 0);
metadata.setPixelsSizeC(toInt(sizeC), 0);
metadata.setPixelsSizeT(toInt(sizeT), 0);

% Set channels ID and samples per pixel
for i = 1: sizeC
    metadata.setChannelID(['Channel:0:' num2str(i-1)], 0, i-1);
    metadata.setChannelSamplesPerPixel(toInt(1), 0, i-1);
end

%%%%%%%%%%%%%%
%ROI
%
% parameters.roi.mask: a 1*n structure that contains n mask files.
%
% parameters.roi.mode: a numeric vector, where 0 - defaut (), 1 - custom by parameters.roi.design_matrix. Indicate how each mask should be processed.
%
% parameters.roi.design_matrix: a numeric matrix, rows are masks, cols are its associates c,t,z,
% -1 means apply to all such dimension, which equals to matlab's index operator ':'.
% For example, first row: [-1,1,2], which means add first mask to c=:,t=1,c=2.
% Note: masks order correspond to masks whose mode 1. Masks that has mode 0 won't
% show up in design matrix.

if isfield( parameters, 'roi' )
    %icaoberg 2020-03-03 the function ml_ls returns on files in the current
    %dir if the input variable is empty
    if isempty( parameters.roi.masks )
        masks = [];
    else
        masks=ml_ls(parameters.roi.masks);
    end
    mode_1_count=0;
    for i=1:size(masks,2)
        roi_id=['ROI:' int2str(i-1)];
        ROIIndex=i-1;
        ROIRefIndex=i-1;
        metadata.setImageROIRef(roi_id, 0, ROIRefIndex); % void setImageROIRef(String roi, int imageIndex, int ROIRefIndex)
        metadata.setROIID(roi_id, ROIIndex); %setROIID(String id, int ROIIndex)
        metadata.setPolygonID('0', ROIIndex, 0); %setPolygonID(String id, int ROIIndex, int shapeIndex)
        %% mask2roi
        img_mask=tif2img(masks{i});
        % img_mask_re=imresize(img_mask,[sizeX sizeY]); %%%% ???? x 2, y 1
        % B = bwboundaries(img_mask_re);
        B = bwboundaries(img_mask');
        ss=size(B{1});
        points_str='';
        for ii=1:ss(1)
            points_str=[points_str num2str(B{1}(ii,1)) ',' num2str(B{1}(ii,2)) ' '];
        end
        %%
        %setPolygonPoints(String points, int ROIIndex, int shapeIndex)
        metadata.setPolygonPoints(java.lang.String(['points0[' points_str(1:end-1) ']']), ROIIndex, 0);
        if parameters.roi.mode(i)==0
            continue
        elseif parameters.roi.mode(i)==1
            mode_1_count=mode_1_count+1;
            roi_C=parameters.roi.design_matrix(mode_1_count,1);
            roi_T=parameters.roi.design_matrix(mode_1_count,2);
            roi_Z=parameters.roi.design_matrix(mode_1_count,3);
            if roi_C~=-1
                %setPolygonTheC(NonNegativeInteger theC, int ROIIndex, int shapeIndex)
                metadata.setPolygonTheC(ome.xml.model.primitives.NonNegativeInteger(java.lang.Integer(roi_C)),ROIIndex,0);
            end
            if roi_T~=-1
                metadata.setPolygonTheT(ome.xml.model.primitives.NonNegativeInteger(java.lang.Integer(roi_T)),ROIIndex,0);
            end
            if roi_Z~=-1
                metadata.setPolygonTheZ(ome.xml.model.primitives.NonNegativeInteger(java.lang.Integer(roi_Z)),ROIIndex,0);
            end
        else
            warning('Wrong parameters.roi.mode. parameters.roi.mode should be 0/1 as int.');
            return;
        end
    end
end
%MapAnnotation
% parameters.MapAnnotation: a n*2 structure that contains n MapAnnotation pairs.

if isfield( parameters, 'MapAnnotation' )
    pairList=java.util.ArrayList;
    for i=1:size(parameters.MapAnnotation,1)
        key=parameters.MapAnnotation{i,1};
        value=parameters.MapAnnotation{i,2};
        pairList.add(ome.xml.model.MapPair(key,value));
    end
    mapAnnotationIndex=0;
    annotationRefIndex=0;
    annotation=['Annotation:' int2str(mapAnnotationIndex)];
    metadata.setMapAnnotationID(annotation,mapAnnotationIndex); %void setMapAnnotationID(String id,int mapAnnotationIndex)
    metadata.setMapAnnotationValue(pairList,mapAnnotationIndex); %void setMapAnnotationValue (List<MapPair> value,int mapAnnotationIndex)
    metadata.setImageAnnotationRef(annotation,0,annotationRefIndex) %setImageAnnotationRef(String annotation, int imageIndex, int annotationRefIndex)
end
%%%%%%%%%%%%%%
end

function dimensionOrders = getDimensionOrders()
% List all values of DimensionOrder
dimensionOrderValues = javaMethod('values', 'ome.xml.model.enums.DimensionOrder');
dimensionOrders = cell(numel(dimensionOrderValues), 1);
for i = 1 :numel(dimensionOrderValues)
    dimensionOrders{i} = char(dimensionOrderValues(i).toString());
end
end

function is = is_octave ()
    is = exist ('OCTAVE_VERSION', 'builtin') == 5;
end