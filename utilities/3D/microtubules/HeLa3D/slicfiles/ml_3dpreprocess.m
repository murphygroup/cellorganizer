function [image,binimage] = ml_3dpreprocess( image, cropimage, ...
			contstr, bgsub, threshmeth)
% [PROCIMAGE,BINIMAGE] = ML_3DPREPROCESS( IMAGE, CROPIMAGE, CONTSTR, ...
%			... BGSUB, THRESHMETH)
% 
% Processes a raw image that will be used in feature extraction. First, 
% hot spots are removed from the image, then the micrograph is contrast 
% stretched, background subtracted, and finally thresholded. 
% 
% Inputs
% 	IMAGE: a raw image of type 'uint8', 'uint16', or 'double'
% 	CROPIMAGE: mask used to isolate single cell patterns
%	CONTSTR: type of contrast stretching of image 
% 		'upper' - stretch max voxel intensity to 255
%		'lower' - shift min voxel intensity to 0
% 		'both' - min-> 0, max-> 255
%		'none' - image is not stretched
%	BGSUB: 	'nobgsub' - do not perform background subtraction
% 		set as anything else to allow background sub.
% 	THRESHMETH: thresholding method
%		'nih' - nih thresholding (default)
%		'rc' - rc thresholding
% Outputs
% 	PROCIMAGE - the processed image,'uint8'
% 	BINIMAGE - the thresholded (binary) image
%
% Adapted by Justin Newberg 2 December 2005

% Copyright (C) 2006  Murphy Lab
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

if ~exist( 'threshmeth','var')
    threshmeth = 'nih';
end

if ~isa(bgsub,'char')
  bgsub = 'yesbgsub';
end

% Remove hot spots from the image
thresh = 0;
CLASS = class(image);

switch CLASS
    case 'uint8'
        [counts bin] = imhist( image(:),2^8);
    case 'uint16'
        [counts bin] = imhist( image(:),2^16);
    case 'double'
        %[counts bin] = hist( image(:),[max(image(:)):1:min(image(:))]);
        [counts bin] = hist( image(:),[min(image(:)):1:max(image(:))]);
        % Yanhua Hu, Feb 07
    otherwise
        disp('error in datatype');
        return;
end

u = find( counts) - 1;
n = 1;
while n<length(u),
    t1 = double(u(n));
    t2 = double(u(n+1));
    if (t2 > 1.5*t1) * (t2 >= 300)
        thresh = t1;
        break;
    end
    n = n+1;
end

if thresh>0
    image(find(image>thresh)) = min(u);
end

% Contrast stretch the image
switch contstr
    case 'both'
        image = double(image);
        imagemin = min(image(:));
        image = uint8(floor((image-imagemin)*255/...
                    (max(image(:))-imagemin)));
    case 'lower'
        image = uint8(floor(double(image - min(image(:)))));
    case 'upper'
        image = uint8(floor(double(image)*255/double(max(image(:)))));
    case 'none'
    otherwise
        error('Unknown contrast stretching method');
end

if ~strcmp( bgsub,'nobgsub')
    image = ml_3dbgsub( image);
end

if(exist('cropimage','var') & (~isempty(cropimage)))
    image = ml_mask(image,cropimage);
end

switch threshmeth
    case 'nih'
        thresh = 255*ml_threshold(image);
    case 'rc'
        thresh = ml_rcthreshold( image);
    otherwise
        error('Unknown thresholding method');
end

binimage = ml_binarize( image, uint8(floor(thresh)));
