function [names, features, slfnames] = ml_3dslfcalc( image, featsetname, ... 
						cropimage, dnaimage, scale, threshmeth)

% ML_3DSLFCALC    Feature calculation on images. Meant to replace 
% 		  ml_3dfeatset.m in NEWSLIC
% [NAMES FEATURES SLFNAMES] = ML_3DSLFCALC( IMAGE, FEATSETNAME, ... 
% 						CROPIMAGE, DNAIMAGE) 
% 
% 	IMAGE = 3D image
% 
% 	FEATSETNAME
% 
% 	      'SLF9'          - SLF9 28 unselected morphological features
% 	      'SLF10'         - SLF10 9 SDA selected features from SLF9
% 	      'SLF11'         - SLF11 42 unselected morphological, edge and 
% 				Haralick texture features at 0.5 microns and
%				64 gray levels
% 	      'SLF14'         - SLF14 14 Subset of SLF9 (all non DNA features)
% 	      'SLF17'         - SLF17 7 first SDA selected morphological, 
%				edge, and Haralick texture features at 0.4 
%				microns and 256 gray levels
% 	      'SLF18'         - SLF18 34 first 34 SDA selected features from 
% 				SLF11 (specific for 3D 3T3 Set)
% 	      'SLF19'         - SLF19 56 unselected morphological, edge, and
% 				Haralick texture features and  the 14 DNA 
% 				features from SLF9
% 	      'SLF20'         - SLF20 52 SDA selected features from SLF19
% 	      'SLF22'         - SLF11 + 14 DNA features
% 
% 	CROPIMAGE = binary 2D/3D mask (optional), which defines the region of 
%		    single cell. If not defined, supply empty matrix []
% 
%	DNAIMAGE = 3D DNA image (optional), if not defined, supply empty 
%		   matrix [] 
%
% 	SCALE = the resolution of the raw image. If not defined, this is 
%		set to 1x1x1 um

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

% Created by Matt Joffee and Justin Newberg on 11-23-05

% Check arguments & use defaults where necessary

if ~exist('threshmeth', 'var')
    threshmeth = 'nih';
end

MaxArgs = 5;
RequiredArgs = 1;
if( nargin < RequiredArgs) % too few arguments
    error(['This function takes at least' sprintf('%i',RequiredArgs) ...
            ' arguments']);
end
%if( nargin > MaxArgs) % too many arguments
if( nargin > MaxArgs+1) % Y.H. Feb 07, for the new parameter 'threshmeth'
    error(['This function takes at most' sprintf('%i',RequiredArgs) ...
            ' arguments']);
end

scale_3DHela = [0.05 0.05 0.2];
scale_3D3T3 = [0.11 0.11 0.5];
scale_3DUCE = [0.0977 0.0977 0.1628];

featsetname = upper(featsetname);

preproc.contstr = 'none';
preproc.bgsub = 'yesbgsub';
%preproc.threshmeth = 'nih';
preproc.threshmeth = threshmeth;  % Y.H. Feb 07, threshmeth is now a param.

morph.minobjsize = [1];
morph.scale = [1 1 1];
morph.scaledist = [1 1 203/48.8];

edge.method = 'normal';

text.tratio = [1 1];
text.method = 'summation';
text.tgray = [];

% Define featureset-specific flags
switch featsetname
case 'SLF9'
    featidx = [1:28];
    if exist( 'scale','var') 
        morph.scaledist = scale/scale(1);
    else
        scale = scale_3DHela;
    end
case 'SLF10'
    featidx = [14,24,13,4,2,1,15,21,22];
    if exist( 'scale','var') 
        morph.scaledist = scale/scale(1);
    else
        scale = scale_3DHela;
    end
case 'SLF11'
    featidx = [1:8,15:20,29:56];
    if exist( 'scale','var') 
        morph.scaledist = scale/scale(1); 
    else
        scale = scale_3D3T3;
    end
    resolution = [0.55 0.55 0.5];
    %text.tratio = scale./resolution;		% 9/25/06
    text.tratio = resolution./scale;		% 10/05/06
    text.tratio(1) = [];
    text.tgray = 64;
    edge.method = 'masked';
    morph.minobjsize = 2;
    preproc.contstr = 'upper';
case 'SLF14'
    featidx = [1:8 15:20];
    if ~exist( 'scale','var') 
        morph.scaledist = scale/scale(1);
    else
        scale = scale_3DHela;
    end
case 'SLF17'
    featidx = [4,30,35,41,42,54,56];
    if exist( 'scale','var') 
        morph.scaledist = scale/scale(1);
    else
        scale = scale_3DHela;
    end
    resolution = [0.4 0.4 0.4];
    %text.tratio = scale./resolution;		% 9/25/06
    text.tratio = resolution./scale;		% 10/05/06
    text.tratio(1) = [];
    text.tgray = 256;
    text.method = 'average';
case 'SLF18'
    featidx = [30,33,37,45,5,42,3,35,38,47,36,43, ...
            39,48,2,55,41,40,51,54,49,50,34,46, ...
            4,52,16,15,32,19,6,31,56,18];
    if exist( 'scale','var') 
        morph.scaledist = scale/scale(1); 
    else
        scale = scale_3D3T3;
    end
    resolution = [0.55 0.55 0.5];
    %text.tratio = scale./resolution;		% 9/25/06
    text.tratio = resolution./scale;		% 10/05/06
    text.tratio(1) = [];
    text.tgray = 64;
    edge.method = 'masked';
    morph.minobjsize = 2;
    preproc.contstr = 'upper';
case 'SLF19'
    featidx = [1:56];
    if exist( 'scale','var') 
        morph.scaledist = scale/scale(1);
    else
        scale = scale_3DUCE;
    end
    resolution = scale_3DUCE;
    %text.tratio = scale./resolution;		% 9/25/06
    text.tratio = resolution./scale;		% 10/05/06
    text.tratio(1) = [];
    if ~strcmp(class(image),'uint8')            % 12/10/06 YH
	preproc.contstr = 'upper';
    end
    
case 'SLF20'
    featidx = [1:7,9:19,21:22,24,25,27:56];
    if exist( 'scale','var')
        morph.scaledist = scale/scale(1);
    else
        scale = scale_3DUCE;
    end
    resolution = scale_3DUCE;
    %text.tratio = scale./resolution;		% 9/25/06
    text.tratio = resolution./scale;		% 10/05/06
    text.tratio(1) = [];
case 'SLF22'
    featidx = [1:56];
    if exist( 'scale','var') 
        morph.scaledist = scale/scale(1); 
    else
        scale = scale_3D3T3;
    end
    resolution = [0.55 0.55 0.5];
    %text.tratio = scale./resolution;		% 9/25/06
    text.tratio = resolution./scale;		% 10/05/06
    text.tratio(1) = [];
    text.tgray = 64;
    edge.method = 'masked';
    morph.minobjsize = 2;
    preproc.contstr = 'upper';
otherwise
    error( 'Unrecognized feature set name');
end

% Validate inputs (image)
temp = size(image);

if( (length(temp) ~= 3) | (min(temp) < 2) | (~(isnumeric(image)) ))
    error('input image must be a 3D numeric matrix ');
end

% Validate inputs (dnaimage/cropimage)

if(~(exist('dnaimage','var')))
    dnaimage = [];
end

if(~(exist('cropimage','var')))
    cropimage = [];
end

% Calculate features
[n,f,s] = ml_3dfeats(image,dnaimage,cropimage,featidx,...
		preproc,morph,edge,text); 

% Output features
features = f;
names = n;
slfnames = s;


