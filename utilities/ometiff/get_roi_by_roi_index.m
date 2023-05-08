function img = get_roi_by_roi_index(img_path, channel_index, roi_index )
% ROI = GET_ROI_by_roi_index( IMAGE_PATH, CHANNEL, ROIID )
%
% Input
% * img_path (path to a valid OME.TIFF file)
% * roi_index (index of a ROI in the image)
% * channel_index (index of a channel in the image)
%
% Output
% * The region of interest (ROI) image corresponding to the ID

% Xin Lu (xlu2@andrew.cmu.edu)
%
% Copyright (C) 2017-2018 Murphy Lab
% Computational Biology Department
% Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.o
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
% 09132017 icaoberg Updated method to dilate the edge of the image returned
% by the ROI call to Bioformats
% 5/5/2023 tedz: Update roi reader to process properly formated
% bioformats metadata rois.
% 5/6/2023 R.F. Murphy give formatted display for roi# and #channels


img = [];

try
    reader = bfGetReader(img_path);
    omeMeta = reader.getMetadataStore();
catch err
    warning('This method requires BioFormats for Matlab');
    getReport( err )
    return
end

% Check channel index argument 

channel_count = get_number_of_channels( img_path );
if channel_count < channel_index || channel_index < 1
    warning('Input argument channel_index is out of range.');
    img = tif2img( img_path );
    return
end

number_of_rois = get_number_of_rois(img_path);
if number_of_rois == 0
    warning('Input image contains no ROIs.');
    img = tif2img( img_path );
    return
end

if number_of_rois < roi_index || roi_index < 1
    warning('Input argument roi_index is out of range.');
    img = tif2img( img_path );
    return 
else

disp(['File has ' num2str(channel_count) ' channels; reading channel ' ...
    num2str(channel_index) ' for ROI ' num2str(roi_index)]);

%     Get all the shapes belonging to the ROI at this index
mask = uint8(zeros(str2num(omeMeta.getPixelsSizeY(0)), str2num(omeMeta.getPixelsSizeX(0)))); 

%   The code considers that a ROI can contain multiple shapes, but for now
%   our images only contain ROIs with one shape
for j = 1:omeMeta.getShapeCount(roi_index-1)
   if strcmp(omeMeta.getShapeType(roi_index-1, j-1), 'Polygon') == 1
         % Xin's original code for getting Polygon ROI points
      Points = omeMeta.getPolygonPoints(roi_index-1, j-1);
      points_ = strsplit(char(Points),'points');
      %format to read if points_ is just 1x1 array
      s = size(points_);
      if s(2) == 1
          points_id = points_{1};
          points_id = strsplit(points_id,' ');
      else
          points_id = points_{2};
          f_ind = strfind(points_id,'[');
          points_id = strsplit(points_id(f_ind+1:end-2), ' ');
      end
      for k = 1:numel(points_id)
          point = strsplit(points_id{k}, ',');
          if s(2) == 1
              x = str2num(point{1}) + 1;
              y = str2num(point{1}) + 1;
          else
              x = str2num(point{1});
              y = str2num(point{2});
          end
          mask(x,y) = 1;
      end
      se = strel( 'disk', 2, 4 );
      mask1 = imdilate( mask, se );
      mask = imfill( mask1, 'hole' );
   else
      % Code added for getting Rectangle ROI points
      RHeight = double(omeMeta.getRectangleHeight(roi_index-1, j-1));
      RWidth = double(omeMeta.getRectangleWidth(roi_index-1, j-1));
      RX1 = double(omeMeta.getRectangleX(roi_index-1, j-1));
      RY1 = double(omeMeta.getRectangleY(roi_index-1, j-1));
      RX2 = RX1+RWidth;
      RY2 = RY1+RHeight;
      c = [RX1 RX2 RX2 RX1 RX1];
      r = [RY1 RY1 RY2 RY2 RY1];
      newmask = roipoly(imread(img_path), c, r);
      mask = mask | newmask;
  end 
end

if channel_count == 1
    % This works for blood cell images
    img = imread(img_path);
    mask = uint8(mask);
    mask = cat(3, mask, mask, mask); 
    img = uint8(img).*mask;
else
    % Xin's code for 2D/3D HeLa ROI
    img1 = OME_loadchannel(img_path, channel_index);
    mask = repmat(mask, 1, 1, size(img1, 3));
    img = uint8(img1).*mask;
end
reader.close();
end%get_roi_by_roi_index