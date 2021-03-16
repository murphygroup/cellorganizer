function [names, features, slfnames] = ml_3dfeatset( image, featsetname, ...
						  cropimage,dnaimage, ...
						  tratio,tgray, ...
						  scale,bgsub,threshmeth)

% ML_3DFEATSET    Feature calculation for raw images
% [NAMES, FEATURES, SLFNAMES] = ML_3DFEATSET( IMAGE, FEATSETNAME,CROPIMAGE, 
%						DNAIMAGE, TRATIO, TGRAY, 
%						SCALE, BGSUB)
%
%	IMAGE = 3D image
%
%
%       FEATSETNAME = 
%
%             'SLF9'          - SLF9 28 unselected morphological features 	       
%             'SLF10'         - SLF10 9 SDA selected features from SLF9	       
%             'SLF11'         - SLF11 42 unselected morphological, edge and Haralick texture features	       
%             'SLF14'         - SLF14 14 Subset of SLF9 (all non DNA features)
%             'SLF17'         - SLF17 7 first 7 SDA selected features from SLF11 (specific for 3D HeLa Set)
%             'SLF18'         - SLF18 34 first 34 SDA selected features from SLF11 (specific for 3D 3T3 Set)	       
%             'SLF19'         - SLF19 56 SLF11 & the 14 DNA features from SLF9
%             'SLF20'         - SLF20 52 SDA selected features from SLF19 
%
%	CROPIMAGE = binary 2D/3D mask (optional), which defines the region of single cell
%		    if not defined supply empty matrix []
%	
%	DNAIMAGE = 3D DNA image (optional), if not defined supply empty matrix []
%
%	TRATIO = A vector specifying 3D downsapling ratio for texture 
%		feature calculation. It should be a scalar or a 2 
%		element vector. All values should be greater than or 
%		equal to 1. If it is a scalar, all 3 dimensions will be 
%		down sampled at this ratio (TRATIO). If a 2 element 
%		vector is provided, the first element stands for the 
%		ratio on the X-Y plane and the second element stands for 
%		the ratio on Z. 
% 
%	TGRAY = Number of gray levels used in texture feature 
%		calculation. 1 <= TGRAY <= 256. UINT8
%
% 	SCALE = micrometers per pixel. Only used if TRATIO is empty & 
%		FEATSETNAME = 'SLF17' or 'SLF18'. If a 2 element vector 
%		is provided, the first element stands for the length and 
%		width in micrometers of the pixel in the X-Y plane and 
%		the second element stands for the depth of the pixel in 
%		the Z plane. Otherwise, scale is a scalar just for the 
%		X-Y plane.
%
% 	BGSUB should be omitted if you want background subtraction, or
% 		'nobgsub' if you do not want it.
%

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

% Created by Matthew Joffe on 5-23-05 
% "Happy 3D Feature calculating" :)
%

% Check arguments & use defaults where necessary
MaxArgs = 9;
RequiredArgs = 2;
if( nargin < RequiredArgs) % too few arguments
    error(['This function takes at least' sprintf('%i',RequiredArgs) ...
            ' arguments']);
end
if( nargin > MaxArgs) % too many arguments
    error(['This function takes at least' sprintf('%i',RequiredArgs) ...
            ' arguments']);
end

if ~exist('threshmeth','var')
    threshmeth='nih';
end

if ~exist('scale','var')
    scale = [];
end

if isempty(scale)
    scale=[1 1 1];
end

% Sort out which features need calculating
switch featsetname
    
case 'SLF9'
    featidx = [1:28];
case 'SLF10'
    featidx = [14,24,13,4,2,1,15,21,22];
case 'SLF11'
    featidx = [1:8,15:20,29:56];
case 'SLF14'
    featidx = [1:8,15:20];
case 'SLF17'
    featidx = [4,30,35,41,42,54,56];
case 'SLF18'
    featidx = [30,33,37,45,5,42,3,35,38,47,36,43, ...
            39,48,2,55,41,40,51,54,49,50,34,46, ...   
            4,52,16,15,32,19,6,31,56,18];
case 'SLF19'
    featidx = [1:56];
case 'SLF20'
    featidx = [1:7,9:19,21:22,24,25,27:56];
otherwise
    error( 'Unrecognized feature set name');
end


% Validate inputs (image)

temp = size(image);

if( (length(temp) ~= 3) | (min(temp) < 2) | (~(isnumeric(image)) ))
    error('input image must be a 3D numeric matrix ');
end

% Validate inputs (dnaimage)

if(~(exist('dnaimage','var')))
    dnaimage = [];
end



% Validate inputs (tratio,scale)

if(~(exist('tratio','var')) | isempty(tratio))
    if( strcmp(featsetname,'SLF17') | strcmp(featsetname,'SLF18') )
        if(exist('scale','var') & (length(scale(:)) <= 2) & isnumeric(scale))
            if(strcmp(featsetname,'SLF17'))
                tratio = repmat(0.4,[1,length(scale)])./scale;
            else
                tratio = repmat(0.5,[1,length(scale)])./scale;
            end
            if(min(tratio)<1)
                error('invalid scale -> tratio conversion : tratio < 1');
            end
        else
            error('scale invalid');
        end
    else
        tratio = [];
    end
else
    if( ~( (length(tratio(:)) <= 2) & (min(tratio) >= 1)) )
        error('invalid tratio');
    end
end

% Validate inputs (tgray)

if(exist('tgray','var')&(~(isempty(tgray))))
    if( ~( (length(tgray) < 2) & (isnumeric(tgray)) & (tgray  <= 256 ) & (tgray > 1.51) ) )
        error('invalid tgray');
    else
        tgray = round(tgray);
    end	
else
    tgray = [];
end


% Preprocess the image (convert uint8)

image = double(image);
imagemax=max(image(:));
imagemin=min(image(:));
image=floor((image-imagemin)*255/(imagemax-imagemin));
image = uint8(image);


% Preprocess the image (background subtraction)

if ~exist('bgsub','var')
    bgsub='yesbgsub';
end

if ~strcmp(bgsub,'nobgsub')
    image = ml_3dbgsub(image);
end


% Preprocess the image (masking)

if(exist('cropimage','var') & (~isempty(cropimage)))
    image=ml_mask(image,cropimage);
end


% Calculate features

[n,f,s] = ml_3dfeat(image,dnaimage,featidx,tratio,tgray,scale,threshmeth);

% Output features

features = f;
names = n;
slfnames = s;


