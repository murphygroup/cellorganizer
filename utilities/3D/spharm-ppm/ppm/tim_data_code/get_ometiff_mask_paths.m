function masks_paths=get_ometiff_mask_paths(imgs_paths,intermediate_output_path)
masks_paths={};
mask_folder=[intermediate_output_path filesep 'masks'];
if ~exist(mask_folder)
    system(['mkdir -p ' mask_folder]);
end
for i=1:length(imgs_paths)
    mask=get_mask_from_ometiff(imgs_paths{i});
    mask_path=[mask_folder filesep 'mask' int2str(i) '.tif'];
    masks_paths{end+1}=mask_path;
    img2tif(mask,mask_path);
end
end

function mask = get_mask_from_ometiff(img_path)

%%%% TO BE EDITED:
% ROI = GET_ROI( IMAGE_PATH, CHANNEL, ROIID )
%
% Input
% * img (a valid OME.TIFF file)
% * roiid (a valid ROI id)
%
% Output
% * The region of interest (ROI) image corresponding to the ID
%%%% END

% Xin Lu (xlu2@andrew.cmu.edu)
%
% Copyright (C) 2017-2019 Murphy Lab
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

img = [];

try
    reader = bfGetReader(img_path);
    omeMeta = reader.getMetadataStore();
catch err
    warning('This method requires BioFormats for Matlab');
    getReport( err )
    return
end

number_of_rois=get_number_of_rois(img_path);
if number_of_rois==0
    warning( ['No ROIs found in ' img_path] );
    img=tif2img( img_path );
    return
end

channel=1;
%retrieve the number of channels
if omeMeta.getChannelCount(0)<channel
    error('Selected channel index out of bounds.');
end

roiid=1;
if number_of_rois<roiid
    error('ROI ID out of ROI number!');
else
    %get ROI points
    Points=omeMeta.getPolygonPoints(0,0);
    points=char(Points);
    points_=strsplit(points,'points');
    points_id=points_{roiid+1};
    f_ind=strfind(points_id,'[');
    points_id=points_id(f_ind+1:end-2);
    points_id=strsplit(points_id,' ');
    
    mask=uint8(zeros(str2num(omeMeta.getPixelsSizeX(0)),str2num(omeMeta.getPixelsSizeY(0))));
    for i=1:numel(points_id)
        point=strsplit(points_id{i},',');
        % mask(str2num(point{1}),str2num(point{2}))=1;
        mask(str2num(point{2}),str2num(point{1}))=1;
    end
    
    se = strel('disk',2,4);
    mask=imdilate( mask, se );
    mask=imfill(mask,'hole');
end

% img = OME_loadchannel(img_path,channel);
% mask=repmat(mask, 1, 1, size(img,3));
% img = uint8(img).*mask;
end