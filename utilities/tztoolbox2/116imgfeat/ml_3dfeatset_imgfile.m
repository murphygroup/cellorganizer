function [names,features,slfnames] = ...
    ml_3dfeatset_imgfile(imgfiles,featsetname,cropfiles,dnafiles,varargin)
%ML_3DFEATSET_IMGFILE Cacluate SLF for image files.
%   [NAMES,FEATURES] = 
%       ML_3DFEATSET_IMGFILE(IMGFILES,FEATSETNAME,CROPFILES,DNAFILES)
%   returns SLF of image files in the cell array IMGFILES, in which
%   each element is the path of an image file. CROPFILES and DNAFILES
%   are files for cropping masks and dna images. CROPFILES and DNAFILES
%   could be empty, which means there are no such image and features
%   will be only calculated on origial images in IMGFILES.
%
%   [NAMES,FEATURES] = 
%       ML_3DFEATSET_IMGFILE(IMGFILES,FEATSETNAME,CROPFILES, ...
%       DNAFILES,TRATIO,TGRAY,SCALE,BGSUB,THRESHMETH)
%   also considers other parameters. 
%
%   [NAMES,FEATURES,SLFNAMES] = ML_3DFEATSET_IMGFILE(...) returns SLF
%   names.
%   
%   For more details for undescribe inputs and outputs, 
%   see ML_3DFEATSET.

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

%   01-Aug-2005 Initial write TINGZ Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 4
    error('At least 4 arguments are required')
end

if isempty(imgfiles)
    error('No images are input');
end

if isempty(featsetname)
    error('Please specify a valid feature set name');
end

for i=1:length(imgfiles)
    img(:,:,i)=ml_readimage(imgfiles{i});
end

if isempty(cropfiles)
    cropimage=[];
else
    if (length(cropfiles)~=1) & (length(cropfiles)~=length(imgfiles))
        error('wrong mask files');
    end
    for i=1:length(cropfiles)
        cropimage(:,:,i)=ml_readimage(cropfiles{i});
    end
    
end

if isempty(dnafiles)
    dnaimage=[];
else
    for i=1:length(dnafiles)
        dnaimage(:,:,i)=ml_readimage(dnafiles{i});
    end  
end

[names, features, slfnames] = ml_3dslfcalc( img, featsetname, ...
    cropimage,dnaimage, ...
    varargin{:});
